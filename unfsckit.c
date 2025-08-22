/* campbell, 2024-2025. isc license */
/* minimum viable usage: printf 'hello\n' | ./fsckit | ./unfsckit */
#include "unfsckit.h"
#include "fft_anywhere.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <string.h>

__attribute((weak)) unsigned char verbose = 1;

static float cmagsquaredf(const float complex x) {
    return crealf(x) * crealf(x) + cimagf(x) * cimagf(x);
}

static float complex renormalize(const float complex x) {
    /* assuming x is already near unity, renormalize to unity w/o div or sqrt */
    return x * (3.0f - cmagsquaredf(x)) * 0.5f;
}

struct argmax {
    float value;
    float power;
};

static struct argmax circular_argmax_of_complex_vector(const size_t S, const float complex s[restrict static S]) {
    float max = 0;
    size_t is_max = 0;
    for (size_t is = 0; is < S; is++) {
        const float this = cmagsquaredf(s[is]);
        if (this > max) {
            max = this;
            is_max = is;
        }
    }

    if (max <= 1.0f / FLT_MAX)
        return (struct argmax) { .value = FLT_MAX };

    /* TODO: use three-point exact dft expr instead of quadratic fit to log of magsquared */
    const float one_over_this = 1.0f / max;
    const float prev = cmagsquaredf(s[(is_max + S - 1) % S]);
    const float next = cmagsquaredf(s[(is_max + 1) % S]);

    const float alpha = logf(prev * one_over_this);
    const float gamma = logf(next * one_over_this);
    const float p = 0.5f * (alpha - gamma) / (alpha + gamma);

    return (struct argmax) {
        .value = remainderf(is_max + p, S),
        .power = max
    };
}

static float complex cosisinf(const float x) {
    return cosf(x) + I * sinf(x);
}

static void dechirp(const size_t S, const size_t L, const size_t H,
                    float complex fft_input[restrict static S],
                    const float window[restrict static S],
                    const float complex history[restrict static H], const size_t ih,
                    const float complex advances[restrict static S * L], const char down,
                    const float freq_offset) {
    /* extract critically sampled values from history, and multiply by conjugate of chirp */
    float complex carrier = 1.0f;

    /* TODO: optimize this */
    const float complex extra = cosisinf(2.0f * (float)M_PI * freq_offset / S);

    for (size_t is = 0; is < S; is++) {
        /* TODO: document this indexing math wow */
        fft_input[is] = history[(is * L + ih) % H] * (down ? carrier : conjf(carrier)) * window[is];
        carrier = renormalize(carrier * advances[(is * L) % (S * L)] * extra);
    }
}

static void populate_advances(const size_t S, const size_t L,
                              float complex advances[restrict static S * L]) {
    /* construct a lookup table of the S*L roots of unity we need for dechirping the input.
     these will be uniformly spaced about the unit circle starting near -1 on the real axis */
    float complex advance = -1.0f * cosisinf(2.0f * (float)M_PI * 0.5f / S);;
    const float complex advance_advance = cosisinf(2.0f * (float)M_PI / (S * L));
    for (size_t isl = 0; isl < S * L; isl++) {
        advances[isl] = advance;
        advance = renormalize(advance * advance_advance);
    }
}

static void butterworth_biquads(const size_t S, float num[restrict static S][3],
                                float den[restrict static S][3],
                                const float fs, const float fc) {
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

static float soft_bit_decision_from_fft(const size_t ibit, const size_t S,
                                        const float complex s[restrict static S]) {
    /* returns +1.0 if the bit is definitely clear, -1.0 if definitely set, 0 if erasure */
    float power_if_clr = 0.0f, power_if_set = 0.0f;
    for (size_t is = 0; is < S; is++)
    /* if this bit is set in the gray code of this symbol... */
        if (gray(is) & (1U << ibit))
            power_if_set += cmagsquaredf(s[is]);
        else
            power_if_clr += cmagsquaredf(s[is]);

    const float denominator = power_if_clr + power_if_set;
    return denominator ? (power_if_clr - power_if_set) / denominator : 0.0f;
}

static unsigned char hamming(unsigned char x) {
    return x | ((((x >> 0U) ^ (x >> 1U) ^ (x >> 2U)) & 0x1) << 4U |
                (((x >> 1U) ^ (x >> 2U) ^ (x >> 3U)) & 0x1) << 5U |
                (((x >> 0U) ^ (x >> 1U) ^ (x >> 3U)) & 0x1) << 6U);
}

struct soft_decode_result {
    unsigned char iword;
    float score;
};

static struct soft_decode_result soft_decode_hamming_naive(const size_t interleave,
                                                    const float soft_bit_history[restrict 7 * interleave]) {
    /* perform naive exhaustive ML decoding of 7 soft bits to get 4 hard bits */
    float best = 0;
    unsigned char iword_best = 0;

    /* loop over all 2^K unique length-N code words */
    for (unsigned char iword = 0; iword < 16; iword++) {
        const unsigned char h = hamming(iword);

        /* take dot product between soft bits and the hamming code of this code word */
        float acc = 0;
        for (unsigned char ibit = 0, mask = 1; ibit < 7; ibit++, mask <<= 1) {
            /* assumes that +1.0 means the bit is clear, -1.0 means the bit is set */
            const float y = soft_bit_history[interleave * ibit];
            if (h & mask) acc -= y;
            else acc += y;
        }

        if (acc > best) {
            best = acc;
            iword_best = iword;
        }
    }

    if (verbose >= 2) dprintf(2, "%s: best score %.2f\r\n", __func__, best);

    return (struct soft_decode_result) {
        .iword = iword_best,
        .score = best
    };
}

static unsigned long long soft_decode_block_hamming(const size_t interleave,
                                                    const float soft_bit_history[restrict 7 * interleave]) {
    /* this takes a block of 7 * interleave soft bits and returns 4 * interleave hard bit
     decisions, giving some robustness against burst interference in the chirp layer */
    /* TODO: replace interleaved hamming codes with something more intelligent */
    float score_sum = 0.0f;
    unsigned long long ret = 0;

    for (size_t iblock = 0; iblock < interleave; iblock++) {
        struct soft_decode_result res = soft_decode_hamming_naive(interleave, soft_bit_history + iblock);
        ret |= res.iword << 4U * iblock;
        score_sum += res.score;
    }
    if (verbose >= 1)
        dprintf(2, "%s: best scores sum to %.2f\r\n", __func__, score_sum);
    return ret;
}

/* these are not really knobs, just magic numbers tied to the specific hamming code */
#define HAMMING_N 7
#define HAMMING_K 4

/* KNOB: half of the number of butterworth filter poles. overall cpu usage is linear with
 this parameter. too many stages can distort the passband response, too few admits noise */
#define BIQUAD_STAGES 2

void unfsckit(const int16_t * (* get_next_sample_func)(const int16_t **, size_t *, void *), void * get_ctx,
              void (* packet_success_function)(const size_t, const unsigned char *, const size_t, void *), void * put_ctx,
              void (* preamble_detected_func)(const size_t, void *), void * preamble_detected_ctx,
              const float sample_rate, const float f_carrier, const float bandwidth,
              const unsigned bits_per_sweep, const unsigned interleave) {
    const size_t S = 1U << bits_per_sweep; /* number of unique measurable symbols */
    /* critically sampled also means there are S complex samples of the band per symbol */

    /* KNOB: intermediate oversampling factor. must be a power of 2, but can be 1.
     using a value larger than 1 allows finer time alignment of input to the demodulator,
     at the expense of more sram usage (but NOT more computational load). the demodulator
     itself always runs at the critical sampling rate s.t the chirps exactly wrap around */
    const size_t L = 4;
    const float sample_rate_filtered = L * bandwidth;

    /* initial value and advance rate of the local oscillator for basebanding */
    float complex carrier = 1.0f;
    const float complex advance = cosisinf(2.0f * (float)M_PI * f_carrier / sample_rate);

    /* compute filter coefficients for butterworth biquad cascade */
    /* KNOB: this scaling factor on bandwidth in the filter parameters strongly affects the
     result. higher values admit more noise, lower values distort the passband */
    float biquad_num[BIQUAD_STAGES][3], biquad_den[BIQUAD_STAGES][3];
    butterworth_biquads(BIQUAD_STAGES, biquad_num, biquad_den, sample_rate, 1.5f * bandwidth);

    const float input_samples_per_filtered_sample = sample_rate / sample_rate_filtered;

    const size_t bytes_expected_max = 256;
    unsigned char * bytes = malloc(bytes_expected_max);

    /* length of buffer of critically sampled complex timeseries */
    const size_t H = 4 * S * L;

    /* properties of the block code */
    const size_t N = HAMMING_N * interleave, K = HAMMING_K * interleave;

    struct planned_forward_fft * plan = plan_forward_fft_of_length(S);
    float complex * restrict const basebanded_ring = malloc(sizeof(float complex) * H);
    float complex * restrict const fft_input = malloc(sizeof(float complex) * S);
    float complex * restrict const fft_output = malloc(sizeof(float complex) * S);
    float complex * restrict const advances = malloc(sizeof(float complex) * S * L);
    float * restrict const soft_bit_history = malloc(sizeof(float) * N);

    /* construct a periodic Hann window */
    float * restrict const window = malloc(sizeof(float) * S);
    for (size_t is = 0; is < S; is++)
        window[is] = 1.0f - cosf(2.0f * (float)M_PI * is / S);

    populate_advances(S, L, advances);

    /* biquad cascade filter state */
    float complex vprev[BIQUAD_STAGES][2] = { 0 };

    /* internal state for simple high pass filtering prior to the biquad cascade */
    int16_t sample_prev = 0;

    /* decimator state */
    float input_samples_since_filtered_sample = 0;

    /* index of what stage of the state machine we are in, so that we can structure this
     as a simple nonthreatening function call from the perspective of calling code */
    unsigned char state = 0;

    /* writer and reader cursors for input ring buffer. the reader shifts its own cursor
     around during preamble detection to align with symbol frames */
    size_t isample_decimated = 0, isample_decimated_next_frame = S * L;

    /* for estimating the absolute timing of things to pass downstream */
    size_t isample_decimated_preamble_start = 0;

    /* after preamble detection, this will hold a residual (circular) shift, in units
     of bins, that should be applied to the fft argmax to get the encoded value */
    float freq_offset = 0;

    /* counters needed by data states */
    unsigned bytes_expected = 0, ibyte = 0;

    /* this will be (re)initialized and then updated as each data symbol comes in */
    unsigned hash = 0;

    unsigned long long decoded_bits = 0;
    unsigned char decoded_bits_filled = 0;

    size_t ih_bit = 0;

    unsigned iframe = 0, iframe_at_last_reset = 0;

    /* buffer of last four upsweep shifts, for preamble detection */
    float prior_upsweeps[4] = { 0, 0, 0, 0 };

    /* this loop structure consists of an outer loop over chunks of samples (of initially
     unknown size), and an inner loop over individual samples within the chunk. this allows
     the same decoder to be used in hard-realtime microcontroller contexts or soft-realtime
     unix pipe type contexts, by replacing the callbacks provided to this function */
    size_t stride;
    for (const int16_t * samples, * end; (samples = get_next_sample_func(&end, &stride, get_ctx));)
        for (; samples != end; samples += stride) {
            const int16_t sample = *samples;
            /* multiply incoming real-valued sample by local oscillator for basebanding */
            float complex filtered = (sample - sample_prev) * conjf(carrier);
            sample_prev = sample;
            carrier = renormalize(carrier * advance);

            /* apply the biquad filter stages to reject stuff outside the passband */
            for (size_t is = 0; is < BIQUAD_STAGES; is++)
                filtered = cfilter(filtered, vprev[is], biquad_num[is], biquad_den[is]);

            /* decimate the filter output to a small integer multiple of nyquist rate */
            input_samples_since_filtered_sample++;
            if (input_samples_since_filtered_sample + 0.5f < input_samples_per_filtered_sample) continue;
            input_samples_since_filtered_sample -= input_samples_per_filtered_sample;

            write(3, &filtered, sizeof(float complex));

            /* store the basebanded filtered decimated samples in a ring buffer */
            basebanded_ring[(isample_decimated++) % H] = filtered;

            /* wait for the buffer to be aligned with the next expected frame */
            if (isample_decimated != isample_decimated_next_frame) continue;
            isample_decimated_next_frame += S * L;

            if (0 == state) {
                /* resetting everything */
                iframe_at_last_reset = iframe;
                freq_offset = 0;
                ih_bit = 0;
                decoded_bits = 0;
                decoded_bits_filled = 0;
                ibyte = 0;
                state++;
            }

            /* if listening for preamble... */
            if (1 == state) {
                /* apply a bunch of criteria for advancing out of preamble detection state */
                do {
                    const float ref = prior_upsweeps[(iframe + 1) % 4];
                    if (FLT_MAX == ref) break;
                    if (FLT_MAX == prior_upsweeps[(iframe + 0) % 4]) break;

                    const float mean_of_middle_upsweeps = (remainderf(prior_upsweeps[(iframe + 0) % 4] - ref, S) + ref +
                                                           remainderf(prior_upsweeps[(iframe + 1) % 4] - ref, S) + ref) * 0.5f;

                    /* not enought time elapsed to see preamble upsweeps */
                    if (iframe - iframe_at_last_reset < 5) break;

                    /* upsweeps do not agree on their value */
                    if (fabsf(remainderf(prior_upsweeps[(iframe + 0) % 4] - ref, S) + ref - mean_of_middle_upsweeps) >= 0.5f ||
                        fabsf(remainderf(prior_upsweeps[(iframe + 1) % 4] - ref, S) + ref - mean_of_middle_upsweeps) >= 0.5f)
                        break;

                    /* assuming freq error is within +/- S/4, the possible time error is
                     bound to +/- S/4 of this value */
                    const float time_offset_midpoint = mean_of_middle_upsweeps;
                    const int shift_midpoint = lrintf(L * time_offset_midpoint);

                    if (verbose >= 2)
                        dprintf(2, "%s: frame %u: two upsweeps agree %.7f %.7f -> %.3f, implied alignment if zero frequency shift is %zu\r\n", __func__, iframe,
                                remainderf(prior_upsweeps[(iframe + 0) % 4] - ref, S) + ref,
                                remainderf(prior_upsweeps[(iframe + 1) % 4] - ref, S) + ref,
                                mean_of_middle_upsweeps, isample_decimated - shift_midpoint);

                    if (verbose >= 2)
                        dprintf(2, "%s: frame %u: decimated sample %zu, considering whether %zu to %zu contains a downsweep\r\n", __func__,
                            iframe, isample_decimated,
                            isample_decimated - S * L - S * L / 2 - shift_midpoint,
                            isample_decimated - S * L - S * L / 2 + S * L - shift_midpoint);

                    dechirp(S, L, H, fft_input, window, basebanded_ring, isample_decimated - S * L - S * L / 2 - shift_midpoint, advances, 0, -0.5f * S - shift_midpoint / (float)L);
                    fft_evaluate_forward(fft_output, fft_input, plan);
                    const struct argmax argmax_up_test = circular_argmax_of_complex_vector(S, fft_output);

                    dechirp(S, L, H, fft_input, window, basebanded_ring, isample_decimated - S * L - S * L / 2 - shift_midpoint, advances, 1, -0.5f * S - shift_midpoint / (float)L);
                    fft_evaluate_forward(fft_output, fft_input, plan);
                    const struct argmax argmax_dn_test = circular_argmax_of_complex_vector(S, fft_output);

                    /* tested interval less likely to be a downsweep than an upsweep */
                    if (argmax_dn_test.power <= argmax_up_test.power) break;

                    if (verbose >= 2)
                        dprintf(2, "%s: frame %u: downsweep detected (%.2f > %.2f)\r\n", __func__, iframe, 10.0f * log10f(argmax_dn_test.power), 10.0f * log10f(argmax_up_test.power));

                    /* got a downsweep, should be able to unambiguously resolve time and
                     frequency shift now, as long as |frequency shift| < bw/4 */
                    const float dn_test_wrapped = -(remainderf(-argmax_dn_test.value - mean_of_middle_upsweeps, S) + mean_of_middle_upsweeps);

                    /* these are both in units of critical samples or frequency bins */
                    const float time_offset = remainderf(0.5f * (mean_of_middle_upsweeps - dn_test_wrapped) - time_offset_midpoint, S) + time_offset_midpoint;

                    /* time offset falls outside of possible range */
                    if (time_offset < time_offset_midpoint - 0.25f * S ||
                        time_offset > time_offset_midpoint + 0.25f * S) {
                        /* TODO: does this ever happen or is this check redundant? */
                        dprintf(2, "%s: frame %u: rejecting due to time alignment violation\r\n", __func__, iframe);
                        break;
                    }

                    freq_offset = remainderf(mean_of_middle_upsweeps - time_offset, S);

                    const float shift_unquantized = L * time_offset;
                    const int shift = lrintf(shift_unquantized);

                    /* TODO: freq_offset should be corrected for time residual error */

                    if (verbose >= 2)
                        dprintf(2, "%s: frame %u: time offset %.3f bins (%.3f samples), freq offset %.3f bins (%.3f Hz)\r\n", __func__, iframe,
                                time_offset, time_offset * L, freq_offset, freq_offset * bandwidth / S);

                    /* consider the prior frame as both an upsweep and downsweep */
                    dechirp(S, L, H, fft_input, window, basebanded_ring, isample_decimated - 2 * S * L - shift, advances, 0, freq_offset);
                    fft_evaluate_forward(fft_output, fft_input, plan);
                    const struct argmax argmax_up_prior = circular_argmax_of_complex_vector(S, fft_output);

                    dechirp(S, L, H, fft_input, window, basebanded_ring, isample_decimated - 2 * S * L - shift, advances, 1, freq_offset);
                    fft_evaluate_forward(fft_output, fft_input, plan);
                    const struct argmax argmax_dn_prior = circular_argmax_of_complex_vector(S, fft_output);

                    if (verbose >= 2)
                        dprintf(2, "%s: frame %u: decimated sample %zu: considering whether %zu to %zu contains a downsweep\r\n", __func__,
                                iframe, isample_decimated, isample_decimated - 2 * S * L - shift, isample_decimated - 2 * S * L - shift + S * L);

                    /* tested interval less likely to be a downsweep than an upsweep */
                    if (argmax_dn_prior.power <= argmax_up_prior.power) break;

                    /* if it was a downsweep with the expected shift... */
                    if (fabsf(remainderf(argmax_dn_test.value + shift_unquantized / (float)L - argmax_dn_prior.value + freq_offset, S)) >= 1.0f) break;

                    if (verbose >= 2)
                        dprintf(2, "%s: frame %u: current and previous frame both downsweeps %.3f %.3f\r\n", __func__,
                                iframe, argmax_dn_test.value + shift_unquantized / (float)L, argmax_dn_prior.value - freq_offset);
                    isample_decimated_next_frame -= shift;

                    /* note: we do not refine the residual using the first of the
                     two downchirps as it is contaminated by filter ringing */

                    if (verbose >= 1)
                        dprintf(2, "%s: frame %u: preamble detected, data frame starts at time %zu, implied carrier offset %.2f Hz\r\n",
                                __func__, iframe, isample_decimated_next_frame - S * L, (freq_offset * bandwidth / S));

                    isample_decimated_preamble_start = (isample_decimated_next_frame - 7 * S * L + L / 2) / L - 1;

                    if (preamble_detected_func)
                    /* estimate the absolute time of the start of the first chirp, in units of 1/bandwidth */
                        preamble_detected_func(isample_decimated_preamble_start, preamble_detected_ctx);

                    state++;
                } while(0);

                /* if we did not advance out of detection state... */
                if (state < 2) {
                    /* retrieve one chirp worth of stuff from the buffer, and de-upsweep it */
                    dechirp(S, L, H, fft_input, window, basebanded_ring, isample_decimated - S * L, advances, 0, 0.0f);

                    /* do an fft of the dechirped symbol frame */
                    fft_evaluate_forward(fft_output, fft_input, plan);

                    /* find index (incl estimating the fractional part) of the loudest fft bin */
                    prior_upsweeps[iframe % 4] = circular_argmax_of_complex_vector(S, fft_output).value;
                }
            } else {
                /* retrieve one chirp worth of stuff from the buffer, and de-upsweep it */
                dechirp(S, L, H, fft_input, window, basebanded_ring, isample_decimated - S * L, advances, 0, freq_offset);

                /* do an fft of the dechirped symbol frame */
                fft_evaluate_forward(fft_output, fft_input, plan);

                /* find index (incl estimating the fractional part) of the loudest fft bin */
                const struct argmax argmax = circular_argmax_of_complex_vector(S, fft_output);

                /* KNOB: scaling factor by which our running estimate of the residual error
                 is nudged by each new symbol. longer symbols want higher values of this
                 knob, as there are fewer update opportunities */
                freq_offset += 0.25f * (argmax.value - lrintf(argmax.value));

                /* TODO: possibly update time alignment according to trend in residual */

                for (size_t ibit = 0; ibit < bits_per_sweep; ibit++) {
                    soft_bit_history[ih_bit++] = soft_bit_decision_from_fft(ibit, S, fft_output);

                    /* whenever we have a complete block of interleaved hamming codes... */
                    if (ih_bit == N) {
                        const unsigned long long bits = soft_decode_block_hamming(interleave, soft_bit_history);

                        decoded_bits |= bits << decoded_bits_filled;
                        decoded_bits_filled += K;

                        while (decoded_bits_filled >= 8) {
                            const unsigned char byte = decoded_bits;
                            decoded_bits >>= 8;
                            decoded_bits_filled -= 8;

                            if (0 == ibyte)
                            /* the first byte encodes the size of the message, from 1 to 256 */
                                bytes_expected = byte + 1;
                            else if (1 == ibyte) {
                                /* the next byte encodes a checksum of the size */
                                const unsigned char len_hash = (2166136261U ^ (bytes_expected - 1)) * 16777619U;
                                if (bytes_expected > bytes_expected_max ||
                                    byte != len_hash) {
                                    dprintf(2, "%s: length failed check or length, resetting\r\n", __func__);
                                    state = 0;
                                } else {
                                    dprintf(2, "%s: reading %u bytes\r\n", __func__, bytes_expected);

                                    /* initial value for fnv-1a checksum */
                                    hash = 2166136261U;
                                }
                            }
                            else if (ibyte < bytes_expected + 2) {
                                /* update fnv-1a hash of data bytes */
                                hash = (hash ^ byte) * 16777619U;

                                bytes[ibyte - 2] = byte;
                            }
                            else if (ibyte == bytes_expected + 2) {
                                /* use low bits of hash as checksum */
                                const unsigned hash_low_bits = hash & 0xff;

                                dprintf(2, "%s: parity received: %u, calculated %u, %s\r\n", __func__,
                                        byte, hash_low_bits, byte == hash_low_bits ? "pass" : "fail");

                                if (byte == hash_low_bits)
                                    packet_success_function(bytes_expected, bytes, isample_decimated_preamble_start, put_ctx);

                                /* reset and wait for next packet */
                                state = 0;
                            }
                            ibyte++;
                        }

                        ih_bit = 0;
                    }
                }
            }
            iframe++;
        }

    destroy_planned_forward_fft(plan);
    free(window);
    free(advances);
    free(fft_output);
    free(fft_input);
    free(basebanded_ring);
    free(bytes);
    free(soft_bit_history);
}

/* everything that follows is an optional cmdline interface for standalone usage */

static const int16_t * get_samples_from_stdin(const int16_t ** end_p, size_t * stride_p, void * ctx) {
    int16_t * buf = ctx;
    const ssize_t ret = fread(buf, sizeof(int16_t), 32, stdin);
    if (ret <= 0) return NULL;
    *end_p = buf + ret;
    *stride_p = 1;
    return buf;
}

static int last_byte = -1;

static void write_payload_to_stdout(const size_t B, const unsigned char * bytes, const size_t it, void * ctx) {
    (void)ctx;
    (void)it;
    fwrite(bytes, 1, B, stdout);
    last_byte = bytes[B - 1];
}

static void preamble_detected(const size_t it, void * ctx) {
    const float bandwidth = *(float *)ctx;
    dprintf(2, "%s: assumed start time %.3f s into timeseries\r\n", __func__, it / bandwidth);
}

#include <unistd.h>
#include <assert.h>

__attribute((weak))
int main(const int argc, char * argv[]) {
    /* get a progname for diagnostic text, regardless of runtime environment */
    const char * o, * progname = argc ? ((o = strrchr(argv[0], '/')) ? o + 1 : argv[0]) : __func__;

    float f_carrier = 1500.0f; /* Hz */
    float bandwidth = 250.0f; /* Hz */
    float sample_rate = 48000.0f; /* samples per second */

    unsigned bits_per_sweep = 5;
    unsigned interleave = 6;

    /* handle --[key]=[value] or space-separated [key] [value] argument pairs */
    for (int iarg = 1; iarg < argc; iarg++) {
        char * key = argv[iarg] + (!strncmp(argv[iarg], "--", 2) ? 2 : 0);
        const size_t keylen = strcspn(key, "=");
        const char * val = '=' == key[keylen] ? key + keylen + 1 : argv[iarg + 1] ? argv[++iarg] : "";
        key[keylen] = '\0';

        if (!strcmp(key, "fs")) sample_rate = strtof(val, NULL);
        else if (!strcmp(key, "fc")) f_carrier = strtof(val, NULL);
        else if (!strcmp(key, "bandwidth")) bandwidth = strtof(val, NULL);
        else if (!strcmp(key, "bits_per_sweep")) bits_per_sweep = strtoul(val, NULL, 10);
        else if (!strcmp(key, "interleave")) interleave = strtoul(val, NULL, 10);

        else {
            dprintf(2, "error: %s: unrecognized argument \"%s\"\n", progname, key);
            exit(EXIT_FAILURE);
        }
    }

    assert(interleave * HAMMING_K <= 56);
    assert(bits_per_sweep > 1);

    setvbuf(stdin, NULL, _IONBF, 0);
    setvbuf(stdout, NULL, _IONBF, 0);

    int16_t buf[32];
    unfsckit(get_samples_from_stdin, buf,
             write_payload_to_stdout, NULL,
             preamble_detected, &(float) { bandwidth },
             sample_rate, f_carrier, bandwidth, bits_per_sweep, interleave);

    if (last_byte != -1 && last_byte != '\n' && isatty(STDOUT_FILENO))
        fputc('\n', stdout);
}
