/* minimum viable usage: printf 'hello\n' | ./fsckit | ./unfsckit */
#include "unfsckit.h"
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

    return (is_max + p) - (is_max + p >= 0.5f * S ? S : 0);
}

static float complex cosisinf(const float x) {
    return cosf(x) + I * sinf(x);
}

static void dechirp(const size_t S, const size_t L,
                    float complex fft_input[restrict static S],
                    const float complex history[restrict static S * L], const unsigned ih,
                    const float complex advances[restrict static S * L], const char down,
                    const float residual) {
    /* extract critically sampled values from history, and multiply by conjugate of chirp */
    float complex carrier = 1.0f;

    /* TODO: optimize this */
    const float complex extra = cosisinf(2.0f * (float)M_PI * residual / S);

    for (size_t is = 0; is < S; is++) {
        /* TODO: document this indexing math wow */
        fft_input[is] = history[(is * L + ih) % (S * L)] * (down ? carrier : conjf(carrier));
        carrier = renormalize(carrier * advances[(is * L) % (S * L)] * extra);
    }
}

static void populate_advances(const size_t S, const size_t L, float complex advances[restrict static S * L]) {
    /* construct a lookup table of the S*L roots of unity we need for dechirping the input.
     these will be uniformly spaced about the unit circle starting near -1 on the real axis */
    float complex advance = -1.0f * cosisinf(2.0f * (float)M_PI * 0.5f / S);;
    const float complex advance_advance = cosisinf(2.0f * (float)M_PI / (S * L));
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
        const float complex apole = wc * cosisinf((float)M_PI * (2.0f * (is + 1) + P - 1.0f) / (2 * P));

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

static unsigned gray(unsigned x) {
    return x ^ (x >> 1U);
}

static float soft_bit_decision_from_fft(const size_t ibit, const size_t S, const float complex s[restrict static S]) {
    float power_zero = 0, power_one = 0;
    for (size_t is = 0; is < S; is++)
    /* if this bit is set in the gray code of this symbol... */
        if (gray(is) & (1U << ibit))
            power_one += cmagsquaredf(s[is]);
        else
            power_zero += cmagsquaredf(s[is]);

    const float denominator = power_zero + power_one;
    return denominator ? (power_zero - power_one) / denominator : 0.0f;
}

static unsigned char hamming(unsigned char x) {
    return x | ((((x >> 0U) ^ (x >> 1U) ^ (x >> 2U)) & 0x1) << 4U |
                (((x >> 1U) ^ (x >> 2U) ^ (x >> 3U)) & 0x1) << 5U |
                (((x >> 0U) ^ (x >> 1U) ^ (x >> 3U)) & 0x1) << 6U);
}

static unsigned char soft_decode_hamming_naive(const size_t ih_bit, const float soft_bit_history[restrict static 32]) {
    /* this could be replaced with well anything */
    float best = 0;
    unsigned char is_best = 0;

    for (unsigned char is = 0; is < 16; is++) {
        const unsigned char h = hamming(is);
        float acc = 0;
        for (unsigned char ib = 0, m = 1; ib < 7; ib++, m <<=1) {
            const float y = soft_bit_history[(ih_bit + ib) % 32];
            if (h & m) acc -= y;
            else acc += y;
        }

        if (acc > best) {
            best = acc;
            is_best = is;
        }
    }
    return is_best;
}

void unfsckit(const int16_t * (* get_next_sample_func)(const int16_t **, size_t *, void *), void * get_ctx,
                     void (* put_bytes_func)(const unsigned char *, const size_t, void *), void * put_ctx,
                     const float sample_rate, const float f_carrier, const float bandwidth,
                     const unsigned bits_per_sweep) {
    const size_t S = 1U << bits_per_sweep; /* number of unique measurable symbols */
    /* critically sampled also means there are S complex samples of the band per symbol */

    /* optional intermediate oversampling factor. must be a power of 2, but can be 1. using
     a value larger than 1 allows finer time alignment of the input to the demodulator, at
     the expense of more sram usage (but not more computational load). the demodulator
     itself always runs at the critical sampling rate s.t the chirps exactly wrap around */
    const size_t L = 2;
    const float sample_rate_filtered = L * bandwidth;

    /* initial value and advance rate of the local oscillator for basebanding */
    float complex carrier = 1.0f;
    const float complex advance = cosisinf(2.0f * (float)M_PI * f_carrier / sample_rate);

    /* compute filter coefficients for four-pole butterworth biquad cascade */
    float num[4][3], den[4][3];
    butterworth_biquads(num, den, 4, sample_rate, 0.75f * bandwidth);

    const float input_samples_per_filtered_sample = sample_rate / sample_rate_filtered;

    const size_t bytes_expected_max = 256;
    unsigned char * bytes = malloc(bytes_expected_max);

    struct planned_forward_fft * plan = plan_forward_fft_of_length(S);
    float complex * restrict const history = malloc(sizeof(float complex) * S * L);
    float complex * restrict const fft_input = malloc(sizeof(float complex) * S);
    float complex * restrict const fft_output = malloc(sizeof(float complex) * S);
    float complex * restrict const advances = malloc(sizeof(float complex) * S * L);
    float * restrict const soft_bit_history = malloc(sizeof(float) * 32);

    populate_advances(S, L, advances);

    /* biquad cascade filter state */
    float complex vprev[4][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 } };

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
    unsigned bytes_expected = 0, ibyte = 0;

    /* this will be (re)initialized and then updated as each data symbol comes in */
    unsigned hash;

    size_t ih_bit = 0, ih_bit_used = 0;

    int16_t sample_prev = 0;

    size_t stride;
    for (const int16_t * samples, * end; (samples = get_next_sample_func(&end, &stride, get_ctx));)
        for (; samples != end; samples += stride) {
            const int16_t sample = *samples;
            /* multiply incoming real-valued sample by local oscillator for basebanding */
            float complex filtered = (sample - sample_prev) * conjf(carrier);
            sample_prev = sample;
            carrier = renormalize(carrier * advance);

            /* apply the biquad filter stages to reject stuff outside the passband */
            for (size_t is = 0; is < 4; is++)
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
            dechirp(S, L, fft_input, history, ih, advances, 0, residual);

            /* do an fft of the dechirped symbol frame */
            fft_evaluate_forward(fft_output, fft_input, plan);

            /* find index (incl estimating the fractional part) of the loudest fft bin */
            float power = 0;
            const float value = circular_argmax_of_complex_vector(&power, S, fft_output);

            if (!power) continue;

            /* if listening for preamble... */
            if (0 == state) {
                /* if three or more agreeing upsweeps have been detected, also listen for downsweeps */
                float power_dn = 0, value_dn = FLT_MAX;
                if (upsweeps >= 3) {
                    dechirp(S, L, fft_input, history, ih, advances, 1, -residual);
                    fft_evaluate_forward(fft_output, fft_input, plan);
                    value_dn = circular_argmax_of_complex_vector(&power_dn, S, fft_output);
                }

                /* if not yet detecting downsweeps, or upsweep was notably louder than possible downsweep... */
                if (power >= 2.0f * power_dn) {
                    /* got an upsweep */
                    downsweeps = 0;

                    const float shift_unquantized = value * L;
                    const int shift = lrintf(shift_unquantized);
                    ih_next_frame -= shift;
                    residual += (shift_unquantized - shift) / L;

                    if (fabsf(value) >= 0.5f)
                        upsweeps = 1;
                    else {
                        upsweeps++;
                        /* TODO: do a running average of the residual over all preamble
                         upsweeps instead of just using the most recent value */

                        if (upsweeps >= 2 && shift) fprintf(stderr, "%s: shifting by %d\r\n", __func__, shift);

                        fprintf(stderr, "%s: upsweep at %u: %ld mB, %.2f, total now %u\r\n", __func__, (unsigned)(ih_next_frame - S * L), lrintf(1e3 * log10f(power)), value, upsweeps);
                    }
                } else {
                    /* got a downsweep */
                    downsweeps++;

                    fprintf(stderr, "%s: downsweep at %u: %ld mB, %.2f, total now %u\r\n", __func__, (unsigned)(ih_next_frame - S * L), lrintf(1e3 * log10f(power_dn)), value_dn, downsweeps);

                    if (1 == downsweeps) {
                        /* the value detected here allows us to disambiguate between
                         timing error and carrier offset error, and correct both */

                        const float shift_unquantized = 0.5f * value_dn * L;
                        const int shift = lrintf(shift_unquantized);

                        ih_next_frame += shift;
                        residual += (float)shift / L;

                        fprintf(stderr, "%s: shifted by %d, next frame at %u, residual is %.2f\r\n", __func__,
                                shift, (unsigned)ih_next_frame, residual);
                    }
                    else if (2 == downsweeps && fabsf(value_dn) < 2.0f)
                        state++;
                    else if (2 == downsweeps) {
                        /* just reset and go back to listening for upsweeps */
                        upsweeps = 0;
                        downsweeps = 0;
                        residual = 0;
                    }
                }
            } else {
                /* nudge residual toward error in this bit */
                residual += 0.5f * (value - lrintf(value));

                for (size_t ibit = 0; ibit < bits_per_sweep; ibit++)
                    soft_bit_history[(ih_bit++) % 32] = soft_bit_decision_from_fft(ibit, S, fft_output);

                if (ih_bit - ih_bit_used >= 14) {
                    const unsigned char lo_bits = soft_decode_hamming_naive(ih_bit_used + 0, soft_bit_history);
                    const unsigned char hi_bits = soft_decode_hamming_naive(ih_bit_used + 7, soft_bit_history);
                    ih_bit_used += 14;
                    const unsigned char byte = lo_bits | hi_bits << 4U;

                    if (1 == state) {
                        bytes_expected = byte + 1;
                        if (bytes_expected > bytes_expected_max)
                            bytes_expected = 0;

                        fprintf(stderr, "%s: reading %u bytes\r\n", __func__, bytes_expected);

                        /* initial value for djb2 checksum */
                        hash = 5381;
                        ibyte = 0;
                        state++;
                    }
                    else if (2 == state) {
                        /* update djb2 hash of data bytes */
                        hash = hash * 33U ^ byte;

                        fprintf(stderr, "%s: %u/%u: %u, err %ld ppt\r\n", __func__, ibyte, bytes_expected, byte, lrintf((value - lrintf(value)) * 1e3f));

                        bytes[ibyte++] = byte;

                        if (bytes_expected == ibyte)
                            state++;
                    }
                    else if (3 == state) {
                        /* use low bits of djb2 hash as checksum */
                        const unsigned hash_low_bits = hash & 0xff;

                        fprintf(stderr, "%s: parity received: %u, calculated %u, %s\r\n", __func__,
                                byte, hash_low_bits, byte == hash_low_bits ? "pass" : "fail");

                        if (byte == hash_low_bits)
                            put_bytes_func(bytes, bytes_expected, put_ctx);

                        /* reset and wait for next packet */
                        state = 0;
                        upsweeps = 0;
                        downsweeps = 0;
                        residual = 0;
                        ih_bit = 0;
                        ih_bit_used = 0;
                    }
                }
            }
        }

    destroy_planned_forward_fft(plan);
    free(advances);
    free(fft_output);
    free(fft_input);
    free(history);
    free(bytes);
    free(soft_bit_history);
}

static const int16_t * get_samples_from_stdin(const int16_t ** end_p, size_t * stride_p, void * ctx) {
    int16_t * buf = ctx;
    const ssize_t ret = fread(buf, sizeof(int16_t), 32, stdin);
    if (ret <= 0) return NULL;
    *end_p = buf + ret;
    *stride_p = 1;
    return buf;
}

static void put_bytes_to_stdout(const unsigned char * bytes, const size_t B, void * ctx) {
    (void)ctx;
    fwrite(bytes, 1, B, stdout);
}

__attribute((weak))
int main(void) {
    const unsigned bits_per_sweep = 5;

    /* input arguments, all in cycles, samples, or symbols per second */
    const float sample_rate = 46875, f_carrier = 3662.11f, bandwidth = 366.211;

    setvbuf(stdin, NULL, _IONBF, 0);

    int16_t buf[32];
    unfsckit(get_samples_from_stdin, buf,
             put_bytes_to_stdout, NULL,
             sample_rate, f_carrier, bandwidth, bits_per_sweep);

    fputc('\n', stdout);
}
