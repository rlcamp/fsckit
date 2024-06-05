#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>
#include <string.h>

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
    float fs = 46875.0f;

    /* bandwidth of modulation, also sample rate if critically sampled */
    float bw = 366.211;

    /* centre frequency */
    float fc = 3662.11f;

    /* this parameter also controls the spreading factor */
    unsigned bits_per_sweep = 5;

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

        unsigned bits = 0;
        unsigned short bits_filled = 0;

        /* djb2 */
        unsigned hash = 5381;

        /* emit four unshifted upsweeps, with continuous carrier phase across sweeps */
        carrier = emit_sweep(carrier, T, advances, 0, 0);
        carrier = emit_sweep(carrier, T, advances, 0, 0);
        carrier = emit_sweep(carrier, T, advances, 0, 0);
        carrier = emit_sweep(carrier, T, advances, 0, 0);
        carrier = emit_sweep(carrier, T, advances, 0, 0);
        carrier = emit_sweep(carrier, T, advances, 0, 0);

        /* two unshifted downsweeps */
        carrier = emit_sweep(carrier, T, advances, 0, 1);
        carrier = emit_sweep(carrier, T, advances, 0, 1);

        /* some upsweeps, circularly shifted to encode length of message in bytes */
        bits = hamming_one_full_byte(B - 1);
        bits_filled = 14;

        /* one shifted upsweep per data symbol */
        for (size_t ibyte = 0, hash_enqueued = 0; bits_filled || !hash_enqueued; ) {
            while (bits_filled >= bits_per_sweep) {
                const unsigned symbol = bits & (S - 1U);
                carrier = emit_symbol(carrier, T, advances, symbol, L);
                bits >>= bits_per_sweep;
                bits_filled -= bits_per_sweep;
            }
            if (hash_enqueued && bits_filled) {
                /* emit an extra symbol to complete the last byte if there are leftover bits */
                const unsigned symbol = bits & ((1U << bits_filled) - 1);
                carrier = emit_symbol(carrier, T, advances, symbol, L);
                bits >>= bits_filled;
                bits_filled = 0;
            }
            if (bits_filled < 14 && ibyte < B) {
                const unsigned char byte = bytes[ibyte++];
                hash = hash * 33U ^ byte;
                bits |= hamming_one_full_byte(byte) << bits_filled;
                bits_filled += 14;
            }
            else if (bits_filled < 14 && ibyte == B && !hash_enqueued) {
                /* one byte for checksum (lowest bits of djb2 of data bytes) */
                bits |= hamming_one_full_byte(hash) << bits_filled;
                bits_filled += 14;
                hash_enqueued = 1;
            }
        }

        bytes += B;
    }

    for (size_t ioffset = 0; ioffset < T; ioffset++)
        fwrite(&(int16_t) { 0 }, sizeof(int16_t), 1, stdout);

    free(advances);
}
