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

static float cmagsquaredf(const float complex x) {
    return crealf(x) * crealf(x) + cimagf(x) * cimagf(x);
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

    const float alpha = logf(prev * one_over_this);
    const float gamma = logf(next * one_over_this);
    const float p = 0.5f * (alpha - gamma) / (alpha + gamma);

    return (is_max + p) - (is_max + p >= 0.5f * S ? S : 0);
}

static float complex cosisinf(const float x) {
    return cosf(x) + I * sinf(x);
}

static void dechirp(const size_t S, const size_t L, const size_t H,
                    float complex fft_input[restrict static S],
                    const float complex history[restrict static H], const unsigned ih,
                    const float complex advances[restrict static S * L], const char down,
                    const float residual) {
    /* extract critically sampled values from history, and multiply by conjugate of chirp */
    float complex carrier = 1.0f;

    /* TODO: optimize this */
    const float complex extra = cosisinf(2.0f * (float)M_PI * residual / S);

    for (size_t is = 0; is < S; is++) {
        /* TODO: document this indexing math wow */
        fft_input[is] = history[(is * L + ih) % H] * (down ? carrier : conjf(carrier));
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

static unsigned char soft_decode_hamming_naive(const float soft_bit_history[restrict], const size_t stride) {
    float best = 0;
    unsigned char is_best = 0;

    /* loop over all code words */
    for (unsigned char is = 0; is < 16; is++) {
        const unsigned char h = hamming(is);

        /* take dot product between observed and hamming code of this code word */
        float acc = 0;
        for (unsigned char ib = 0, m = 1; ib < 7; ib++, m <<= 1) {
            /* assumes that +1.0 means the bit is clear, -1.0 means the bit is set */
            const float y = soft_bit_history[stride * ib];
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

static unsigned long long soft_decode_block_hamming(const float soft_bit_history[restrict], const size_t interleave) {
    /* TODO: replace interleaved hamming codes with something more intelligent */
    unsigned long long ret = 0;

    for (size_t iblock = 0; iblock < interleave; iblock++)
        ret |= soft_decode_hamming_naive(soft_bit_history + iblock, interleave) << 4U * iblock;
    return ret;
}

/* these are not really knobs, just magic numbers tied to the specific hamming code */
#define HAMMING_N 7
#define HAMMING_K 4

/* KNOB: half of the number of butterworth filter poles. overall cpu usage is linear with
 this parameter. too many stages can distort the passband response, too few admits noise */
#define BIQUAD_STAGES 2

void unfsckit(const int16_t * (* get_next_sample_func)(const int16_t **, size_t *, void *), void * get_ctx,
              void (* put_bytes_func)(const unsigned char *, const size_t, void *), void * put_ctx,
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
    butterworth_biquads(BIQUAD_STAGES, biquad_num, biquad_den, sample_rate, 0.8f * bandwidth);

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
    unsigned ih = 0, ih_next_frame = S * L;

    /* after preamble detection, this will hold a residual (circular) shift, in units
     of bins, that should be applied to the fft argmax to get the encoded value */
    float residual = 0;

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

            /* store the basebanded filtered decimated samples in a ring buffer */
            basebanded_ring[(ih++) % H] = filtered;

            /* wait for the buffer to be aligned with the next expected frame */
            if (ih != ih_next_frame) continue;
            ih_next_frame += S * L;

            /* retrieve one chirp worth of stuff from the buffer, and de-upsweep it */
            dechirp(S, L, H, fft_input, basebanded_ring, ih - S * L, advances, 0, residual);

            /* do an fft of the dechirped symbol frame */
            fft_evaluate_forward(fft_output, fft_input, plan);

            /* find index (incl estimating the fractional part) of the loudest fft bin */
            float power = 0;
            const float value = circular_argmax_of_complex_vector(&power, S, fft_output);

            if (!power) continue;

            if (0 == state) {
                /* resetting everything */
                iframe_at_last_reset = iframe;
                residual = 0;
                ih_bit = 0;
                decoded_bits = 0;
                decoded_bits_filled = 0;
                ibyte = 0;
                state++;
            }

            /* if listening for preamble... */
            if (1 == state) {
                const float mean_of_oldest_upsweeps = (prior_upsweeps[(iframe + 0) % 4] +
                                                       prior_upsweeps[(iframe + 1) % 4] +
                                                       prior_upsweeps[(iframe + 2) % 4]) * (1.0f / 3.0f);

                if (iframe - iframe_at_last_reset >= 5 &&
                    fabsf(prior_upsweeps[(iframe + 0) % 4] - mean_of_oldest_upsweeps) < 0.5f &&
                    fabsf(prior_upsweeps[(iframe + 1) % 4] - mean_of_oldest_upsweeps) < 0.5f &&
                    fabsf(prior_upsweeps[(iframe + 2) % 4] - mean_of_oldest_upsweeps) < 0.5f) {
                    dprintf(2, "%s: frame %u: oldest three upsweeps agree %.3f\r\n", __func__, iframe, mean_of_oldest_upsweeps);

                    float power_dn = 0;
                    dechirp(S, L, H, fft_input, basebanded_ring, ih - S * L, advances, 1, 0);
                    fft_evaluate_forward(fft_output, fft_input, plan);
                    const float value_dn_now = circular_argmax_of_complex_vector(&power_dn, S, fft_output);
                    if (power_dn > power) {
                        /* got a downsweep, should be able to unambiguously resolve time and
                         frequency shift now, as long as |frequency shift| < bw/4 */
                        const float shift_unquantized = L * (mean_of_oldest_upsweeps - value_dn_now) * 0.5f;
                        const int shift = lrintf(shift_unquantized);

                        /* best guess of residual so far, not yet counting first downsweep */
                        residual = (3.0f * (mean_of_oldest_upsweeps - shift / (float)L) +
                                    value_dn_now + shift / (float)L) * (1.0f / 4.0f);

                        /* consider the prior frame as both an upsweep and downsweep */
                        float power_up_prior = 0, power_dn_prior = 0;
                        dechirp(S, L, H, fft_input, basebanded_ring, ih - 2 * S * L - shift, advances, 0, residual);
                        fft_evaluate_forward(fft_output, fft_input, plan);
                        circular_argmax_of_complex_vector(&power_up_prior, S, fft_output);

                        dechirp(S, L, H, fft_input, basebanded_ring, ih - 2 * S * L - shift, advances, 1, residual);
                        fft_evaluate_forward(fft_output, fft_input, plan);
                        const float value_dn_prior = circular_argmax_of_complex_vector(&power_dn_prior, S, fft_output) - residual;

                        /* if it was a downsweep with the expected shift... */
                        if (power_dn_prior > power_up_prior && fabsf(value_dn_prior - value_dn_now - shift / (float)L) < 0.5f) {
                            dprintf(2, "%s: frame %u: current and previous frame both downsweeps %.3f %.3f\r\n", __func__,
                                    iframe, value_dn_prior, value_dn_now + shift / (float)L);
                            ih_next_frame -= shift;

                            /* final estimate of carrier offset considers last three upsweeps
                             and both downsweeps equally */
                            residual = (3.0f * (mean_of_oldest_upsweeps - shift / (float)L) +
                                        value_dn_now + shift / (float)L +
                                        value_dn_prior) * (1.0f / 5.0f);

                            dprintf(2, "%s: frame %u: data frame starts at time %u, implied carrier offset %ld Hz\r\n",
                                    __func__, iframe, ih_next_frame, lrintf(residual * bandwidth / S));

                            state++;
                        }
                        else dprintf(2, "%s: frame %u: possible downsweep\r\n", __func__, iframe);
                    }
                }

                prior_upsweeps[iframe % 4] = value;
            } else {
                /* KNOB: scaling factor by which our running estimate of the residual error
                 is nudged by each new symbol. longer symbols want higher values of this
                 knob, as there are fewer update opportunities */
                residual += 0.25f * (value - lrintf(value));

                /* if the residual dechirp alignment error has grown large enough to imply
                 a time offset of more than one sample, shift the time alignment */
                const int shift = lrintf(residual * L);
                if (shift) {
                    ih_next_frame -= shift;
                    residual -= (float)shift / L;
                }

                for (size_t ibit = 0; ibit < bits_per_sweep; ibit++) {
                    soft_bit_history[ih_bit++] = soft_bit_decision_from_fft(ibit, S, fft_output);

                    /* whenever we have a complete block of interleaved hamming codes... */
                    if (ih_bit == N) {
                        const unsigned long long bits = soft_decode_block_hamming(soft_bit_history, interleave);

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
                                    put_bytes_func(bytes, bytes_expected, put_ctx);

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

static void put_bytes_to_stdout(const unsigned char * bytes, const size_t B, void * ctx) {
    (void)ctx;
    fwrite(bytes, 1, B, stdout);
    last_byte = bytes[B - 1];
}

#include <unistd.h>
#include <assert.h>

__attribute((weak))
int main(const int argc, const char * const * const argv) {

    /* input arguments, all in cycles, samples, or symbols per second */
    const float f_carrier = argc > 1 ? strtof(argv[1], NULL) : 1500.0f;
    const float bandwidth = argc > 2 ? strtof(argv[2], NULL) : 250.0f;
    const float sample_rate = argc > 3 ? strtof(argv[3], NULL) : 48000.0f;
    const unsigned bits_per_sweep = argc > 4 ? strtoul(argv[4], NULL, 10) : 5;
    const unsigned interleave = argc > 5 ? strtoul(argv[5], NULL, 10) : 6;

    assert(interleave * HAMMING_K <= 56);

    setvbuf(stdin, NULL, _IONBF, 0);
    setvbuf(stdout, NULL, _IONBF, 0);

    int16_t buf[32];
    unfsckit(get_samples_from_stdin, buf,
             put_bytes_to_stdout, NULL,
             sample_rate, f_carrier, bandwidth, bits_per_sweep, interleave);

    if (last_byte != -1 && last_byte != '\n' && isatty(STDOUT_FILENO))
        fputc('\n', stdout);
}
