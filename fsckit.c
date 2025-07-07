/* campbell, 2024-2025, isc license */
#include "fsckit.h"

#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <limits.h>

static uint64_t xorshift64star(void) {
    /* marsaglia et al. generates 64 bits at a time, the most significant bits are the most random */
    static uint64_t x = 1; /* must be nonzero */
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    return x * 0x2545F4914F6CDD1DULL;
}

static float frand_minus_frand(void) {
    /* generate 64 random bits, of which we will use the most significant 46, in two groups of 23 */
    const uint64_t bits = xorshift64star();

    /* generate two random numbers each uniformly distributed on [1.0f, 2.0f) */
    const union { uint32_t u; float f; } x = { .u = 0x3F800000U | ((bits >> 41) & 0x7FFFFFU) };
    const union { uint32_t u; float f; } y = { .u = 0x3F800000U | ((bits >> 18) & 0x7FFFFFU) };

    /* and subtract them, yielding a triangular distribution on [-1.0f, +1.0f] */
    return x.f - y.f;
}

static float complex renormalize(const float complex x) {
    /* assuming x is already near unity, renormalize to unity w/o div or sqrt */
    const float magsquared = crealf(x) * crealf(x) + cimagf(x) * cimagf(x);
    return x * (3.0f - magsquared) * 0.5f;
}

static float complex cosisinf(const float x) {
    return cosf(x) + I * sinf(x);
}

static float complex emit_sweep(void (* emit_sample_func)(void *, const int16_t),
                                void * emit_sample_ctx,
                                float complex carrier, const size_t T,
                                const float complex modulation_advances[restrict static T],
                                const float complex carrier_advance,
                                const size_t shift, const int down, const float amplitude) {
    float complex modulation = 1.0f;

    for (size_t it = 0; it < T; it++) {
        emit_sample_func(emit_sample_ctx, lrintf(cimagf(carrier * modulation) * amplitude + frand_minus_frand()));
        modulation = renormalize(modulation * modulation_advances[((down ? T - it : it) + shift) % T]);
        carrier = renormalize(carrier * carrier_advance);
    }
    return carrier * modulation;
}

static unsigned degray(unsigned x) {
    for (unsigned m = x; m; m >>= 1, x ^= m);
    return x;
}

static unsigned char hamming(unsigned char x) {
    return x | ((((x >> 0U) ^ (x >> 1U) ^ (x >> 2U)) & 0x1) << 4U |
                (((x >> 1U) ^ (x >> 2U) ^ (x >> 3U)) & 0x1) << 5U |
                (((x >> 0U) ^ (x >> 1U) ^ (x >> 3U)) & 0x1) << 6U);
}

/* these are not really knobs, just magic numbers tied to the specific hamming code */
#define HAMMING_N 7
#define HAMMING_K 4

static unsigned long long hamming_interleaved(const unsigned long long data_bits, const size_t interleave) {
    unsigned long long ret = 0;

    for (size_t ib = 0; ib < interleave; ib++) {
        const unsigned int bits_now = hamming((data_bits >> (ib * HAMMING_K)) & 0xF);

        for (size_t ia = 0; ia < HAMMING_N; ia++) {
            const unsigned long long mask = 1ULL << (ib + interleave * ia);
            if (bits_now & (1ULL << ia))
                ret |= mask;
            else
                ret &= ~mask;
        }
    }
    return ret;
}

float complex fsckit(void (* emit_sample_func)(void *, const int16_t), void * emit_sample_ctx,
                     const float amplitude, const float fs, const float fc, const float bw,
                     const unsigned bits_per_sweep, const unsigned interleave, float complex carrier,
                     const size_t V, const struct iovec iovs[restrict static V]) {
    /* number of unique measurable symbols is 2^bits_per_sweep */
    const size_t S = 1U << bits_per_sweep;

    /* oversampling factor, must be an integer for now */
    const size_t L = lrintf(fs / bw);

    /* output samples per symbol */
    const size_t T = S * L;

    /* sweep rate in Hz per second */
    const float df_dt = bw * fs / T;

    /* properties of the block code */
    const size_t N = HAMMING_N * interleave, K = HAMMING_K * interleave;

    float complex * restrict const modulation_advances = malloc(sizeof(float complex) * T);
    const float complex advance_advance = cosisinf(2.0f * (float)M_PI * df_dt / (fs * fs));
    modulation_advances[0] = cosisinf(2.0f * (float)M_PI * -0.5f * bw / fs);

    /* apply recurrence relation to generate lookup table of carrier advancements */
    for (size_t it = 1; it < T; it++)
        modulation_advances[it] = renormalize(modulation_advances[it - 1] * advance_advance);

    /* initial state of carrier */
    const float complex carrier_advance = cosisinf(2.0f * (float)M_PI * fc / fs);

    /* fnv-1a initial value */
    unsigned hash = 2166136261U;

    /* emit four unshifted upsweeps, with continuous carrier phase across sweeps */
    carrier = emit_sweep(emit_sample_func, emit_sample_ctx, carrier, T, modulation_advances, carrier_advance, 0, 0, amplitude);
    carrier = emit_sweep(emit_sample_func, emit_sample_ctx, carrier, T, modulation_advances, carrier_advance, 0, 0, amplitude);
    carrier = emit_sweep(emit_sample_func, emit_sample_ctx, carrier, T, modulation_advances, carrier_advance, 0, 0, amplitude);
    carrier = emit_sweep(emit_sample_func, emit_sample_ctx, carrier, T, modulation_advances, carrier_advance, 0, 0, amplitude);

    /* two unshifted downsweeps */
    carrier = emit_sweep(emit_sample_func, emit_sample_ctx, carrier, T, modulation_advances, carrier_advance, 0, 1, amplitude);
    carrier = emit_sweep(emit_sample_func, emit_sample_ctx, carrier, T, modulation_advances, carrier_advance, 0, 1, amplitude);

    size_t B = 0;
    for (size_t iv = 0; iv < V; iv++)
        B += iovs[iv].iov_len;

    const unsigned char len_hash = (2166136261U ^ (unsigned char)(B - 1)) * 16777619U;

    unsigned long long data_bits = 0;
    size_t data_bits_filled = 0;

    unsigned long long coded_bits = 0;
    size_t coded_bits_filled = 0;

    const struct iovec * iov = iovs;
    size_t iov_byte = 0;

    for (size_t ibyte = 0; ibyte < B + 3 || coded_bits_filled || data_bits_filled; ) {
        while (coded_bits_filled >= bits_per_sweep) {
            const unsigned symbol = coded_bits & (S - 1U);
            carrier = emit_sweep(emit_sample_func, emit_sample_ctx, carrier, T, modulation_advances, carrier_advance, degray(symbol) * L, 0, amplitude);
            coded_bits >>= bits_per_sweep;
            coded_bits_filled -= bits_per_sweep;
        }

        if (ibyte == B + 3 && coded_bits_filled && !data_bits_filled)
        /* if no more transposed bits are coming and we have a partially completed
         sweep, then just enqueue the difference */
            coded_bits_filled = bits_per_sweep;

        /* if we have enough data bits to do another FEC block... */
        while (coded_bits_filled + N <= sizeof(coded_bits) * CHAR_BIT && data_bits_filled >= K) {
            coded_bits |= hamming_interleaved(data_bits, interleave) << coded_bits_filled;
            data_bits >>= K;
            data_bits_filled -= K;
            coded_bits_filled += N;
        }

        /* if no more data bits and we need to complete an FEC block input, then just
         enqueue the difference in zero bits */
        if (ibyte == B + 3 && data_bits_filled && coded_bits_filled < N)
            data_bits_filled = K;

        /* if we can enqueue another full byte... */
        while (data_bits_filled + 8 <= sizeof(data_bits) * CHAR_BIT && ibyte < B + 3) {
            /* the next byte to send is either the message length, its hash, a data
             byte, or the hash of all the data bytes */
            const unsigned char byte = (0 == ibyte ? B - 1 :
                                        1 == ibyte ? len_hash :
                                        ibyte < B + 2 ? ((unsigned char *)iov->iov_base)[iov_byte++] :
                                        ibyte == B + 2 ? hash : 0);

            if (iov_byte == iov->iov_len) {
                iov++;
                iov_byte = 0;
            }

            data_bits |= (unsigned long long)byte << data_bits_filled;
            data_bits_filled += 8;

            if (ibyte >= 2) hash = (hash ^ byte) * 16777619U;
            ibyte++;
        }
    }

    free(modulation_advances);

    const float initial_imag = cimagf(carrier);
    while (cimagf(carrier) * initial_imag > 0.0f) {
        emit_sample_func(emit_sample_ctx, lrintf(cimagf(carrier) * amplitude + frand_minus_frand()));
        carrier = renormalize(carrier * carrier_advance);
    }

    return carrier;
}

#include <stdio.h>

static void fwrite_sample(void * ctx, const int16_t sample) {
    fwrite(&sample, sizeof(int16_t), 1, ctx);
}

__attribute((weak))
int main(const int argc, const char * const * const argv) {
    /* input arguments, all in cycles, samples, or symbols per second */

    /* centre frequency */
    const float fc = argc > 1 ? strtof(argv[1], NULL) : 1500.0f;

    /* bandwidth of modulation, also sample rate if critically sampled */
    const float bw = argc > 2 ? strtof(argv[2], NULL) : 250.0f;

    /* sample rate */
    const float fs  = argc > 3 ? strtof(argv[3], NULL) : 48000.0f;

    /* this parameter also controls the spreading factor */
    const unsigned bits_per_sweep = argc > 4 ? strtoul(argv[4], NULL, 10) : 5;

    const unsigned interleave = argc > 5 ? strtoul(argv[5], NULL, 10) : 6;

    const float amplitude = argc > 6 ? strtof(argv[6], NULL) : 32766.0f;

    if (interleave * HAMMING_N + bits_per_sweep > sizeof(unsigned long long) * CHAR_BIT) {
        fprintf(stderr, "error: %s: interleave too long\n", __func__);
        exit(EXIT_FAILURE);
    }

    const size_t chirp_period_in_samples = (1U << bits_per_sweep) * lrintf(fs / bw);

    /* enforce oversampling factor is an integer */
    if ((size_t)lrintf((1U << bits_per_sweep) * fs / bw) != chirp_period_in_samples) {
        fprintf(stderr, "%s: sample rate not an exact multiple of bandwidth\n", __func__);
        exit(EXIT_FAILURE);
    }

    /* emit one sweep period of quiet samples */
    for (size_t it = 0; it < chirp_period_in_samples; it++)
        fwrite(&(int16_t) { lrintf(frand_minus_frand()) }, sizeof(int16_t), 1, stdout);

    float complex carrier = 1.0f;

    /* loop over lines of text on stdin, or blocks of 256 bytes */
    char bytes[257];

    while (1) {
        if (!fgets(bytes, sizeof(bytes), stdin)) break;
        const size_t B = strlen(bytes);

        carrier = fsckit(fwrite_sample, stdout, amplitude, fs, fc, bw, bits_per_sweep, interleave, carrier,
                         1, &(struct iovec) { .iov_base = bytes, .iov_len = B });
    }

    /* emit some quiet samples to flush the decoder if being piped directly into it */
    for (size_t it = 0; it < chirp_period_in_samples; it++)
        fwrite(&(int16_t) { lrintf(frand_minus_frand()) }, sizeof(int16_t), 1, stdout);
}
