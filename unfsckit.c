/* minimum viable usage: ./fsckit | ./unfsckit */
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

static void dechirp(const size_t S, const size_t L,
                    float complex fft_input[restrict static S],
                    const float complex history[restrict static S * L], const unsigned ih,
                    const float complex advances[restrict static S * L], const char down) {
    /* extract critically sampled values from history, and multiply by conjugate of chirp */
    float complex carrier = 1.0f;

    for (size_t is = 0; is < S; is++) {
        /* TODO: document this indexing math wow */
        fft_input[is] = history[(is * L + ih) % (S * L)] * (down ? carrier : conjf(carrier));
        carrier = renormalize(carrier * advances[(is * L) % (S * L)]);
    }
}

static void populate_advances(const size_t S, const size_t L, float complex advances[restrict static S * L]) {
    /* construct a lookup table of the S*L roots of unity we need for dechirping the input.
     these will be uniformly spaced about the unit circle starting near -1 on the real axis */
    float complex advance = -1.0f * cexpf(I * 2.0f * (float)M_PI * 0.5f / S);;
    const float complex advance_advance = cexpf(I * 2.0f * (float)M_PI / (S * L));
    for (size_t isl = 0; isl < S * L; isl++) {
        advances[isl] = advance;
        advance = renormalize(advance * advance_advance);
    }
}

static void butterworth_biquads(float num[][3], float den[][3], size_t S, float fs, float fc) {
    /* number of poles must be even */
    const size_t P = 2 * S;

    /* prewarp corner frequency for bilinear transform */
    const float wc = 2.0f * tanf((float)M_PI * fc / fs);

    /* each stage implements a conjugate pair of analog poles */
    for (size_t is = 0; is < P / 2; is++) {
        /* analog butterworth pole. the two poles for this stage are this and its conjugate */
        const float complex apole = wc * cexpf(I * (float)M_PI * (2.0f * (is + 1) + P - 1.0f) / (2 * P));

        /* each analog pole results in one digital pole and one digital zero at -1.0 */
        const float complex dpole = (2.0f - apole) / (2.0f + apole);

        /* polynomial coefficients for pair of digital zeros at -1 */
        num[is][0] = 1.0f;
        num[is][1] = 2.0f;
        num[is][2] = 1.0f;

        /* polynomial coefficients for conjugate pair of digital poles */
        den[is][0] = dpole * conjf(dpole);
        den[is][1] = -2.0f * crealf(dpole);
        den[is][2] = 1.0f;

        /* normalize the set of coefficients for unit gain */
        const float den_sum = den[is][0] + den[is][1] + den[is][2];
        const float den_scale = 1.0f / den[is][0], num_scale = den_scale * den_sum / 4.0f;
        for (size_t ik = 0; ik < 3; ik++) num[is][ik] *= num_scale;
        for (size_t ik = 0; ik < 3; ik++) den[is][ik] *= den_scale;
    }
}

static float complex cfilter(const float complex x, float complex vprev[2],
                             const float num[3], const float den[3]) {
    /* operate on complex input and output with real filter coefficients */
    const float complex v =          x - den[1] * vprev[0] - den[2] * vprev[1];
    const float complex y = num[0] * v + num[1] * vprev[0] + num[2] * vprev[1];

    vprev[1] = vprev[0];
    vprev[0] = v;

    return y;
}

int main(void) {
    const unsigned bits_per_sweep = 5;
    const size_t S = 1U << bits_per_sweep; /* number of unique measurable symbols */
    /* critically sampled also means there are S complex samples of the band per symbol */

    /* input arguments, all in cycles, samples, or symbols per second */
    const float sample_rate = 8000, f_carrier = 2000, bandwidth = 250;

    /* optional intermediate oversampling factor. must be a power of 2, but can be 1. using
     a value larger than 1 allows finer time alignment of the input to the demodulator, at
     the expense of more sram usage (but not more computational load). the demodulator
     itself always runs at the critical sampling rate s.t the chirps exactly wrap around */
    const size_t L = 2;
    const float sample_rate_filtered = L * bandwidth;

    /* initial value and advance rate of the local oscillator for basebanding */
    float complex carrier = 1.0f;
    const float complex advance = cexpf(I * 2.0f * (float)M_PI * f_carrier / sample_rate);

    /* compute filter coefficients for four-pole butterworth biquad cascade */
    float num[2][3], den[2][3];
    butterworth_biquads(num, den, 2, sample_rate, 0.75f * bandwidth);

    const float input_samples_per_filtered_sample = sample_rate / sample_rate_filtered;

    struct planned_forward_fft * plan = plan_forward_fft_of_length(S);
    float complex * restrict const history = malloc(sizeof(float complex) * S * L);
    float complex * restrict const fft_input = malloc(sizeof(float complex) * S);
    float complex * restrict const fft_output = malloc(sizeof(float complex) * S);
    float complex * restrict const advances = malloc(sizeof(float complex) * S * L);

    populate_advances(S, L, advances);

    /* biquad cascade filter state */
    float complex vprev[2][2] = { { 0, 0 }, { 0, 0 } };

    /* decimator state */
    float input_samples_since_filtered_sample = 0;

    /* index of what stage of the state machine we are in, so that we can structure this
     as a simple nonthreatening function call from the perspective of calling code */
    unsigned char state = 0;

    /* writer and reader cursors for input ring buffer. the reader shifts its own cursor
     around during preamble detection to align with symbol frames */
    unsigned ih = 0, ih_next_frame = S * L;

    /* after preamble detection, this will hold a residual (circular) shift, in units
     of bins, that should be applied to the fft argmax to get the encoded value */
    float residual = 0;

    /* counters needed by preamble state */
    unsigned upsweeps = 0, downsweeps = 0;

    /* counters needed by data states */
    unsigned data_symbols_expected = 0, idata = 0;

    /* this will be (re)initialized and then updated as each data symbol comes in */
    unsigned hash;

    float downsweep_prev = 0;

    int16_t sample_prev = 0, sample;
    while (fread(&sample, sizeof(int16_t), 1, stdin)) {
        /* multiply incoming real-valued sample by local oscillator for basebanding */
        float complex filtered = (sample - sample_prev) * conjf(carrier);
        sample_prev = sample;
        carrier = renormalize(carrier * advance);

        /* apply the biquad filter stages to reject stuff outside the passband */
        for (size_t is = 0; is < 2; is++)
            filtered = cfilter(filtered, vprev[is], num[is], den[is]);

        /* decimate the filter output to a small integer multiple of nyquist rate */
        input_samples_since_filtered_sample++;
        if (input_samples_since_filtered_sample + 0.5f < input_samples_per_filtered_sample) continue;
        input_samples_since_filtered_sample -= input_samples_per_filtered_sample;

        /* store the basebanded filtered decimated samples in a ring buffer */
        history[(ih++) % (S * L)] = filtered;

         /* wait for the buffer to be aligned with the next expected frame */
        if (ih != ih_next_frame) continue;
        ih_next_frame += S * L;

        /* retrieve one chirp worth of stuff from the buffer, and de-upsweep it */
        dechirp(S, L, fft_input, history, ih, advances, 0);

        /* do an fft of the dechirped symbol frame */
        fft_evaluate_forward(fft_output, fft_input, plan);

        /* find index (incl estimating the fractional part) of the loudest fft bin */
        float power = 0;
        const float value = remainderf(circular_argmax_of_complex_vector(&power, S, fft_output) - residual, S);

        if (!power) continue;

        /* if listening for preamble... */
        if (0 == state) {
            /* if three or more agreeing upsweeps have been detected, also listen for downsweeps */
            float power_dn = 0, value_dn = FLT_MAX;
            if (upsweeps >= 3) {
                dechirp(S, L, fft_input, history, ih, advances, 1);
                fft_evaluate_forward(fft_output, fft_input, plan);
                value_dn = remainderf(circular_argmax_of_complex_vector(&power_dn, S, fft_output) + residual, S);
            }

            /* if not yet detecting downsweeps, or upsweep was notably louder than possible downsweep... */
            if (power >= 2.0f * power_dn) {
                /* got an upsweep */
                downsweeps = 0;

                const float shift_unquantized = value * L;
                const int shift = lrintf(shift_unquantized);
                if (shift) fprintf(stderr, "%s: shifting by %d\n", __func__, shift);
                ih_next_frame -= shift;
                residual += (shift_unquantized - shift) / L;

                if (fabsf(value) >= 0.5f)
                    upsweeps = 1;
                else
                    upsweeps++;
                    /* TODO: do a running average of the residual over all preamble
                     upsweeps instead of just using the most recent value */

                fprintf(stderr, "%s: upsweep detected at %g, total now %u\n",
                        __func__, value, upsweeps);
            } else {
                /* got a downsweep */
                downsweeps++;
                fprintf(stderr, "%s: downsweep detected at %g + %g = %g\n",__func__,
                        value_dn - residual, residual, value_dn);

                if (2 == downsweeps && fabsf(remainderf(value_dn - downsweep_prev, S)) < 2.0f) {
                    /* the value detected here allows us to disambiguate between
                     timing error and carrier offset error, and correct both */

                    const float shift_unquantized = 0.5f * value_dn * L;
                    const int shift = lrintf(shift_unquantized);

                    ih_next_frame += shift;
                    residual += (float)shift / L;

                    fprintf(stderr, "%s: shifted by %d, residual is %g\n", __func__,
                            shift, residual);

                    state++;
                }
                else if (2 == downsweeps) {
                    /* just reset and go back to listening for upsweeps */
                    upsweeps = 0;
                    downsweeps = 0;
                }
                else downsweep_prev = value_dn;
            }
        } else {
            const unsigned symbol = (lrintf(value + S)) % S;
            fprintf(stderr, "%s: %.2f - %.2f = %.2f -> %u\n", __func__, value + residual, residual, value, symbol);

            if (1 == state) {
                data_symbols_expected = symbol + 1;
                fprintf(stderr, "%s: reading %u values\n", __func__, data_symbols_expected);

                /* initial value for djb2 checksum */
                hash = 5381;
                idata = 0;
                state++;
            }
            else if (2 == state) {
                /* update djb2 hash of data symbols */
                hash = hash * 33U ^ symbol;

                fprintf(stderr, "%s: data symbol %u/%u: %u\n", __func__, idata, data_symbols_expected, symbol);

                idata++;
                if (data_symbols_expected == idata) state++;
            }
            else if (3 == state) {
                /* use low bits of djb2 hash as checksum */
                const unsigned hash_low_bits = hash & (S - 1U);

                fprintf(stderr, "%s: parity received: %u, calculated %u, %s\n", __func__,
                        symbol, hash_low_bits, symbol == hash_low_bits ? "pass" : "fail");

                /* reset and wait for next packet */
                state = 0;
                upsweeps = 0;
                downsweeps = 0;
            }
        }
    }

    destroy_planned_forward_fft(plan);
    free(advances);
    free(fft_output);
    free(fft_input);
    free(history);
}
