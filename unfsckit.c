#include "fft_anywhere.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <string.h>

static float cmagsquaredf(const float complex x) {
    return crealf(x) * crealf(x) + cimagf(x) * cimagf(x);
}

static float complex renormalize(const float complex x) {
    /* assuming x is already near unity, renormalize to unity w/o div or sqrt */
    return x * (3.0f - cmagsquaredf(x)) * 0.5f;
}
static float argmax_of_fft_of_dechirped(float * power_max_p,
                                        const size_t S, const size_t L,
                                        float complex fft_output[restrict static S],
                                        float complex fft_input[restrict static S],
                                        const struct planned_forward_fft * plan,
                                        const float complex history[restrict static S * L],
                                        const size_t ih,
                                        const float complex advances[restrict static S * L],
                                        const int icarrier, const int down) {
    float complex carrier = 1.0f;

    for (size_t is = 0; is < S; is++) {
        /* TODO: document this indexing math wow */
        fft_input[is] = history[(is * L + ih) % (S * L)] * (down ? carrier : conjf(carrier));
        carrier = renormalize(carrier * advances[(is * L + icarrier + S * L + 1) % (S * L)]);
    }

    fft_evaluate_forward(fft_output, fft_input, plan);

    float max = 0;
    size_t iw_max = 0;
    for (size_t iw = 0; iw < S; iw++) {
        const float this = cmagsquaredf(fft_output[iw]);
        if (this > max) {
            max = this;
            iw_max = iw;
        }
    }

    if (!max) return FLT_MAX;

    /* TODO: use exact phase expressions instead of quadratic fit in magnitude squared */
    const float one_over_this = 1.0f / max;
    const float prev = cmagsquaredf(fft_output[(iw_max + S - 1) % S]);
    const float next = cmagsquaredf(fft_output[(iw_max + 1) % S]);

    /* these logarithms are very expensive */
    const float alpha = logf(prev * one_over_this), gamma = logf(next * one_over_this);
    const float p = 0.5f * (alpha - gamma) / (alpha + gamma);

    if (power_max_p) *power_max_p = max;
    return iw_max + p;
}

static int wait_for_frame(size_t * ih_p, size_t * ih_next_frame_p, const size_t S,
                          const size_t L,
                          float complex history[restrict static S * L]) {
    while (*ih_p != *ih_next_frame_p) {
        float complex input;
        if (fread(&input, sizeof(float complex), 1, stdin) <= 0) return -1;
        history[((*ih_p)++) % (S * L)] = input;
    }
    *ih_next_frame_p += S * L;
    return 0;
}

int main(void) {
    const unsigned bits_per_sweep = 5;
    const size_t S = 1U << bits_per_sweep; /* number of unique measurable symbols */
    /* critically sampled also means there are S complex samples of the band per symbol */

    const size_t L = 2;
    struct planned_forward_fft * plan = plan_forward_fft_of_length(S);
    float complex * restrict const history = malloc(sizeof(float complex) * S * L);
    float complex * restrict const fft_input = malloc(sizeof(float complex) * S);
    float complex * restrict const fft_output = malloc(sizeof(float complex) * S);
    float complex * restrict const advances = malloc(sizeof(float complex) * S * L);

    /* assume critically sampled s.t. advance exactly wraps around */
    float complex advance = -1.0f;
    const float complex advance_advance = cexpf(I * 2.0f * (float)M_PI / (S * L));

    for (size_t isl = 0; isl < S * L; isl++) {
        advances[isl] = advance;
        advance = renormalize(advance * advance_advance);
    }

    int icarrier = 0;
    size_t ih = 0, ih_next_frame = S * L;

    /* outermost loop repeats once per packet */
    while (1) {
        /* after preamble detection, this will hold a residual (circular) shift, in units
         of bins, that should be applied to the fft argmax to get the encoded value */
        float residual = 0;

        {
            char eof = 0;
            size_t upsweeps = 0, downsweeps = 0;

            /* loop until expected sequence of upsweeps and downsweeps has been seen */
            while (1) {
                if (-1 == wait_for_frame(&ih, &ih_next_frame, S, L, history)) { eof = 1; break; }

                float power_up = 0, power_dn = 0;
                /* always listen for upsweeps in this state */
                const float value_up = argmax_of_fft_of_dechirped(&power_up, S, L, fft_output, fft_input, plan, history, ih, advances, icarrier, 0);

                /* if three or more upsweeps have been detected, also listen for downsweeps */
                const float value_dn = upsweeps >= 3 ? (argmax_of_fft_of_dechirped(&power_dn, S, L, fft_output, fft_input, plan, history, ih, advances, icarrier, 1)) : FLT_MAX;

                fprintf(stderr, "%s: %g %g\n", __func__, power_up, power_dn);

                /* if we got neither, just keep trying */
                if (value_up >= FLT_MAX && value_dn >= FLT_MAX) continue;

                if (power_up >= 2.0f * power_dn) {
                    const float value_up_wrapped = value_up >= S / 2 ? value_up - S : value_up;
                    fprintf(stderr, "%s: upsweep detected at %g\n", __func__, value_up_wrapped);
                    /* got an upsweep */
                    downsweeps = 0;

                    if (fabsf(value_up_wrapped) >= 0.5f)
                        upsweeps = 1;
                    else {
                        upsweeps++;
                        /* TODO: do a running average of the residual over all preamble
                         upsweeps instead of just using the most recent value */
                        residual = value_up_wrapped;
                    }

                    fprintf(stderr, "%s: %zu upsweeps detected\n", __func__, upsweeps);

                    if (fabsf(value_up_wrapped) >= 0.5f / L) {
                        /* got an upsweep that does not agree with anything previously */
                        fprintf(stderr, "%s: shifting by %ld\n", __func__, lrintf(value_up_wrapped * L));
                        ih_next_frame -= lrintf(value_up_wrapped * L);
                    }
                } else if (upsweeps >= 3) {
                    const float value_dn_wrapped = value_dn >= S / 2 ? value_dn - S : value_dn;
                    fprintf(stderr, "%s(%d): got here %g\n", __func__, __LINE__, value_dn_wrapped);
                    /* got a downsweep */
                    if (fabsf(value_dn_wrapped + residual) < 0.5f * S) {
                        downsweeps++;
                        fprintf(stderr, "%s: downsweep detected at %g + %g = %g\n", __func__, value_dn_wrapped, residual, value_dn_wrapped + residual);

                        if (2 == downsweeps) {
                            /* the value detected here allows us to disambiguate between
                             timing error and carrier offset error, and correct both */

                            const float shift_unquantized = 0.5f * (value_dn_wrapped + residual) * L;
                            const int shift = lrintf(shift_unquantized);

                            fprintf(stderr, "%s: shifting by %d\n", __func__, shift);

                            ih_next_frame += shift;
                            residual += shift_unquantized / L;

                            fprintf(stderr, "%s: residual is %g\n", __func__, residual);

                            break;
                        }
                    }
                }
            }

            if (eof) break;
        }

        /* next symbol encodes the frame length */
        if (-1 == wait_for_frame(&ih, &ih_next_frame, S, L, history)) break;

        const float value_header = argmax_of_fft_of_dechirped(NULL, S, L, fft_output, fft_input, plan, history, ih, advances, icarrier, 0) - residual;
        if (value_header >= FLT_MAX) break;
        const size_t data_symbols_expected = (lrintf(value_header + 2 * S) % S) + 1;

        fprintf(stderr, "%s: reading %.2f + 1 -> %zu values\n", __func__, value_header, data_symbols_expected);

        unsigned parity_calculated = 0;

        /* remaining symbols encode the data */
        for (size_t idata = 0; idata < data_symbols_expected; idata++) {
            if (-1 == wait_for_frame(&ih, &ih_next_frame, S, L, history)) break;

            const float value = argmax_of_fft_of_dechirped(NULL, S, L, fft_output, fft_input, plan, history, ih, advances, icarrier, 0) - residual;
            if (value >= FLT_MAX) break;

            const unsigned symbol = lrintf(value + 2 * S) % S;

            /* TODO: checksum that isn't absolutely worthless */
            parity_calculated ^= symbol;

            fprintf(stderr, "%s: data symbol %zu/%zu: %.2f - %.2f = %.2f -> %u\n", __func__, idata, data_symbols_expected, value + residual, residual, value, symbol);

            /* TODO: continue to track the timing offset */
        }

        const float value_parity = argmax_of_fft_of_dechirped(NULL, S, L, fft_output, fft_input, plan, history, ih, advances, icarrier, 0) - residual;
        if (value_parity >= FLT_MAX) break;
        const unsigned parity_received = lrintf(value_parity + S) % S;

        fprintf(stderr, "%s: parity received %u, calculated %u, %s\n", __func__,
                parity_received, parity_calculated, parity_received == parity_calculated ? "pass" : "fail");
    }

    destroy_planned_forward_fft(plan);
    free(advances);
    free(fft_output);
    free(fft_input);
    free(history);
}
