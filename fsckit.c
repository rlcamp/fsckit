#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

static float complex renormalize(const float complex x) {
    /* assuming x is already near unity, renormalize to unity w/o div or sqrt */
    const float magsquared = crealf(x) * crealf(x) + cimagf(x) * cimagf(x);
    return x * (3.0f - magsquared) * 0.5f;
}

static float complex emit_sweep(float complex carrier, const size_t T,
                                const float complex advances[restrict static T],
                                const size_t shift, const int down) {
    for (size_t it = 0; it < T; it++) {
        fwrite(&(int16_t) { lrintf(crealf(carrier) * 32767.0f * 0.25f) }, sizeof(int16_t), 1, stdout);
        carrier = renormalize(carrier * advances[((down ? T - it : it) + shift) % T]);
    }
    return carrier;
}

static unsigned degray(unsigned x) {
    for (unsigned m = x; m; m >>= 1, x ^= m);
    return x;
}

static float complex emit_symbol(float complex carrier, const size_t T,
                        const float complex advances[restrict static T],
                        const unsigned symbol,  const size_t L) {
    return emit_sweep(carrier, T, advances, degray(symbol) * L, 0);
}

static unsigned char hamming(unsigned char x) {
    return x | ((((x >> 0U) ^ (x >> 1U) ^ (x >> 2U)) & 0x1) << 4U |
                (((x >> 1U) ^ (x >> 2U) ^ (x >> 3U)) & 0x1) << 5U |
                (((x >> 0U) ^ (x >> 1U) ^ (x >> 3U)) & 0x1) << 6U);
}

static unsigned hamming_one_full_byte(unsigned char x) {
    return hamming(x & 0xF) | (hamming(x >> 4U) << 7U);
}

int main(void) {
    /* sample rate */
    float fs = 31250;

    /* bandwidth of modulation, also sample rate if critically sampled */
    float bw = 244.141;

    /* centre frequency */
    float fc = 1953.12f;

    /* this parameter also controls the spreading factor */
    unsigned bits_per_sweep = 5;

    const unsigned interleave = 6;

    if (interleave * 7 + bits_per_sweep > sizeof(unsigned long long) * CHAR_BIT) {
        fprintf(stderr, "error: %s: interleave too long\n", __func__);
        exit(EXIT_FAILURE);
    }

    /* number of unique measurable symbols is 2^bits_per_sweep */
    const size_t S = 1U << bits_per_sweep;

    /* oversampling factor, must be an integer for now */
    const size_t L = lrintf(fs / bw);

    /* output samples per symbol */
    const size_t T = S * L;

    /* enforce oversampling factor is an integer */
    assert((size_t)lrintf(S * fs / bw) == T);

    /* sweep rate in Hz per second */
    const float df_dt = bw * fs / T;

    float complex * restrict const advances = malloc(sizeof(float complex) * T);
    const float complex advance_advance = cexpf(I * 2.0f * (float)M_PI * df_dt / (fs * fs));
    advances[0] = cexpf(I * 2.0f * (float)M_PI * (fc - 0.5f * bw) / fs);

    /* apply recurrence relation to generate lookup table of carrier advancements */
    for (size_t it = 1; it < T; it++)
        advances[it] = renormalize(advances[it - 1] * advance_advance);

    /* initial state of carrier */
    float complex carrier = 1.0f;

    /* maybe emit some quiet samples */
    for (size_t ioffset = 0; ioffset < 999; ioffset++)
        fwrite(&(int16_t) { 0 }, sizeof(int16_t), 1, stdout);

    /* loop over lines of text on stdin, emitting one message per line */
    char * bytes = NULL;
    char * line = NULL;
    size_t linecap = 0;
    while (1) {
        if (!bytes || !bytes[0]) {
            if (getline(&line, &linecap, stdin) <= 0) break;
            bytes = line;
        }
        const size_t line_remaining = strlen(bytes);
        const size_t B = line_remaining <= 256 ? line_remaining : 256;

        /* fnv-1a initial value */
        unsigned hash = 2166136261U;

        /* emit four unshifted upsweeps, with continuous carrier phase across sweeps */
        carrier = emit_sweep(carrier, T, advances, 0, 0);
        carrier = emit_sweep(carrier, T, advances, 0, 0);
        carrier = emit_sweep(carrier, T, advances, 0, 0);
        carrier = emit_sweep(carrier, T, advances, 0, 0);

        /* two unshifted downsweeps */
        carrier = emit_sweep(carrier, T, advances, 0, 1);
        carrier = emit_sweep(carrier, T, advances, 0, 1);

        const unsigned char len_hash = (2166136261U ^ (unsigned char)(B - 1)) * 16777619U;

        unsigned long long bits = 0;
        size_t bits_filled = 0;

        unsigned long long bits_transposed = 0;
        size_t bits_transposed_filled = 0;

        const size_t total_hamming_coded_bits = (7 * (B + 3) * 8) / 4;
        const size_t total_interleaved_bits = ((total_hamming_coded_bits + interleave * 7 - 1) / (interleave * 7)) * interleave * 7;
        const size_t total_sweeps = (total_interleaved_bits + bits_per_sweep - 1) / bits_per_sweep;

        for (size_t isweep = 0, ibyte = 0; isweep < total_sweeps; ) {
            while (bits_transposed_filled >= bits_per_sweep) {
                const unsigned symbol = bits_transposed & (S - 1U);
                carrier = emit_symbol(carrier, T, advances, symbol, L);
                bits_transposed >>= bits_per_sweep;
                bits_transposed_filled -= bits_per_sweep;
                isweep++;
            }

            /* TODO: clean up the criteria for flushing this */
            if (isweep + 1 == total_sweeps && bits_transposed_filled)
                /* if no more transposed bits are coming and we have a partially completed
                 sweep, then just enqueue the difference */
                bits_transposed_filled = bits_per_sweep;

            while (bits_transposed_filled + interleave * 7 <= sizeof(bits_transposed) * CHAR_BIT &&
                bits_filled >= interleave * 7) {

                for (size_t ib = 0, ibit = 0; ib < interleave; ib++)
                    for (size_t ia = 0; ia < 7; ia++, ibit++) {
                        const unsigned long long mask = 1ULL << (ib + interleave * ia + bits_transposed_filled);
                        if (bits & (1ULL << ibit))
                            bits_transposed |= mask;
                        else
                            bits_transposed &= ~mask;
                    }

                bits >>= interleave * 7;
                bits_filled -= interleave * 7;
                bits_transposed_filled += interleave * 7;
            }

            /* if we can enqueue another full byte... */
            while (bits_filled + 14 <= sizeof(bits) * CHAR_BIT && ibyte < B + 3) {
                /* the next byte to send is either the message length, its hash, a data
                 byte, or the hash of all the data bytes */
                const unsigned char byte = (0 == ibyte ? B - 1 :
                                            1 == ibyte ? len_hash :
                                            ibyte < B + 2 ? bytes[ibyte - 2] :
                                            ibyte == B + 2 ? hash : 0);
                bits |= (unsigned long long)hamming_one_full_byte(byte) << bits_filled;
                bits_filled += 14;

                if (ibyte >= 2) hash = (hash ^ byte) * 16777619U;

                ibyte++;
            }

            /* if no more bits are coming from upstream and we need to complete a transpose
             group, then just enqueue the difference in zero bits */
            if (ibyte == B + 3 && bits_filled && bits_transposed_filled < interleave * 7)
                bits_filled = interleave * 7;
        }

        bytes += B;
    }

    for (size_t ioffset = 0; ioffset < 4 * T; ioffset++)
        fwrite(&(int16_t) { 0 }, sizeof(int16_t), 1, stdout);

    free(advances);
}
