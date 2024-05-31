/* minimum viable usage: ./fsckit | ./baseband | ./unfsckit */
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

static float logf_good_enough_approx(float x) {
    /* thanks to stef o'rear for pointing out that ieee floats are already a first-order
     approximation of logf when interpreted as an integer, you just have to scale and shift it */
    const union { float f; uint32_t i; } y = { .f = x };
    return y.i * (0.6931472f / 8388608.0f) - 127.0f * 0.6931472f;
}

static float complex renormalize(const float complex x) {
    /* assuming x is already near unity, renormalize to unity w/o div or sqrt */
    return x * (3.0f - cmagsquaredf(x)) * 0.5f;
}

static float circular_argmax_of_complex_vector(float * max_magsquared_p, const size_t S,
                                               const float complex s[restrict static S]) {
    float max = 0;
    size_t is_max = 0;
    for (size_t is = 0; is < S; is++) {
        const float this = cmagsquaredf(s[is]);
        if (this > max) {
            max = this;
            is_max = is;
        }
    }

    if (max_magsquared_p) *max_magsquared_p = max;
    if (!max) return FLT_MAX;

    /* TODO: use three-point exact dft expr instead of quadratic fit to log of magsquared */
    const float one_over_this = 1.0f / max;
    const float prev = cmagsquaredf(s[(is_max + S - 1) % S]);
    const float next = cmagsquaredf(s[(is_max + 1) % S]);

    /* actual logf is expensive, use a fast approx that is good enough for the use case */
    const float alpha = logf_good_enough_approx(prev * one_over_this);
    const float gamma = logf_good_enough_approx(next * one_over_this);
    const float p = 0.5f * (alpha - gamma) / (alpha + gamma);

    return is_max + p;
}

static float argmax_of_fft_of_dechirped(float * power_max_p,
                                        const size_t S, const size_t L,
                                        float complex fft_output[restrict static S],
                                        float complex fft_input[restrict static S],
                                        const struct planned_forward_fft * plan,
                                        const float complex history[restrict static S * L],
                                        const size_t ih,
                                        const float complex advances[restrict static S * L],
                                        const int down) {
    /* extract critically sampled values from history, and multiply by conjugate of chirp */
    float complex carrier = 1.0f;

    for (size_t is = 0; is < S; is++) {
        /* TODO: document this indexing math wow */
        fft_input[is] = history[(is * L + ih) % (S * L)] * (down ? carrier : conjf(carrier));
        carrier = renormalize(carrier * advances[(is * L + S * L + 1) % (S * L)]);
    }

    fft_evaluate_forward(fft_output, fft_input, plan);

    return circular_argmax_of_complex_vector(power_max_p, S, fft_output);
}

static int wait_for_frame(size_t * ih_p, size_t * ih_next_frame_p, const size_t S,
                          const size_t L, float complex history[restrict static S * L]) {
    /* this function ingests new samples, filling a ring buffer long enough to hold one
     sweep at the intermediate (not necessarily critically sampled) sample rate, and returns
     when the buffer is aligned with the next expected frame */
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

    /* optional intermediate oversampling factor. must be a power of 2, but can be 1. using
     a value larger than 1 allows finer time alignment of the input to the demodulator, at
     the expense of more sram usage (but not more computational load). the demodulator
     itself always runs at the critical sampling rate s.t the chirps exactly wrap around */
    const size_t L = 2;
    struct planned_forward_fft * plan = plan_forward_fft_of_length(S);
    float complex * restrict const history = malloc(sizeof(float complex) * S * L);
    float complex * restrict const fft_input = malloc(sizeof(float complex) * S);
    float complex * restrict const fft_output = malloc(sizeof(float complex) * S);
    float complex * restrict const advances = malloc(sizeof(float complex) * S * L);

    /* construct a lookup table of the S*L roots of unity we need for dechirping the input.
     these will be uniformly spaced about the unit circle starting at -1 on the real axis */
    float complex advance = -1.0f;
    const float complex advance_advance = cexpf(I * 2.0f * (float)M_PI / (S * L));
    for (size_t isl = 0; isl < S * L; isl++) {
        advances[isl] = advance;
        advance = renormalize(advance * advance_advance);
    }

    /* writer and reader cursors for input ring buffer. the reader shifts its own cursor
     around during preamble detection to align with symbol frames */
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
                const float value_up = remainderf(argmax_of_fft_of_dechirped(&power_up, S, L, fft_output, fft_input, plan, history, ih, advances, 0) - residual, S);

                /* if two or more agreeing upsweeps have been detected, also listen for downsweeps */
                const float value_dn = upsweeps >= 2 ? remainderf(argmax_of_fft_of_dechirped(&power_dn, S, L, fft_output, fft_input, plan, history, ih, advances, 1) + residual, S) : FLT_MAX;

                /* if we got neither, just keep trying */
                if (!power_up && !power_dn) continue;

                if (power_up >= 2.0f * power_dn) {
                    /* got an upsweep */
                    downsweeps = 0;

                    fprintf(stderr, "%s: upsweep detected at %g, total now %zu\n",
                            __func__, value_up, upsweeps);

                    const float shift_unquantized = value_up * L;
                    const int shift = lrintf(shift_unquantized);
                    if (shift) fprintf(stderr, "%s: shifting by %d\n", __func__, shift);
                    ih_next_frame -= shift;
                    residual += (shift_unquantized - shift) / L;

                    if (fabsf(value_up) >= 0.5f)
                        upsweeps = 1;
                    else {
                        upsweeps++;
                        /* TODO: do a running average of the residual over all preamble
                         upsweeps instead of just using the most recent value */
                    }
                } else if (upsweeps >= 2) {
                    /* got a downsweep */
                    if (fabsf(value_dn) < 0.5f * S) {
                        downsweeps++;
                        fprintf(stderr, "%s: downsweep detected at %g + %g = %g\n",__func__,
                                value_dn - residual, residual, value_dn);

                        if (2 == downsweeps) {
                            /* the value detected here allows us to disambiguate between
                             timing error and carrier offset error, and correct both */

                            const float shift_unquantized = 0.5f * value_dn * L;
                            const int shift = lrintf(shift_unquantized);

                            ih_next_frame += shift;
                            residual += (float)shift / L;

                            fprintf(stderr, "%s: shifted by %d, residual is %g\n", __func__,
                                    shift, residual);

                            break;
                        }
                    }
                }
            }

            if (eof) break;
        }

        unsigned char state = 1;
        size_t data_symbols_expected = 0, idata = 0;

        /* initial value for djb2 checksum */
        unsigned hash = 5381;

        while (1) {
            if (-1 == wait_for_frame(&ih, &ih_next_frame, S, L, history)) break;

            float power = 0;
            const float value = remainderf(argmax_of_fft_of_dechirped(&power, S, L, fft_output, fft_input, plan, history, ih, advances, 0) - residual, S);
            if (!power) break;

            const unsigned symbol = (lrintf(value + S)) % S;
            fprintf(stderr, "%s: %.2f - %.2f = %.2f -> %u\n", __func__, value + residual, residual, value, symbol);

            if (1 == state) {
                data_symbols_expected = symbol + 1;
                fprintf(stderr, "%s: reading %.2f + 1 -> %zu values\n", __func__, value, data_symbols_expected);
                state++;
            }
            else if (2 == state) {
                /* update djb2 hash of data symbols */
                hash = hash * 33U + symbol;

                fprintf(stderr, "%s: data symbol %zu/%zu: %u\n", __func__, idata, data_symbols_expected, symbol);

                idata++;
                if (data_symbols_expected == idata) state++;
            }
            else if (3 == state) {
                const unsigned parity_received = symbol;

                /* use low bits of djb2 hash as checksum */
                const unsigned parity_calculated = hash & (S - 1U);

                fprintf(stderr, "%s: parity received: %u, calculated %u, %s\n", __func__,
                        parity_received, parity_calculated, parity_received == parity_calculated ? "pass" : "fail");
                break;
            }
        }
    }

    destroy_planned_forward_fft(plan);
    free(advances);
    free(fft_output);
    free(fft_input);
    free(history);
}
