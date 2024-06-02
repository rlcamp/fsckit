#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

static float complex renormalize(const float complex x) {
    /* assuming x is already near unity, renormalize to unity w/o div or sqrt */
    const float magsquared = crealf(x) * crealf(x) + cimagf(x) * cimagf(x);
    return x * (3.0f - magsquared) * 0.5f;
}

static float complex emit_sweep(float complex carrier, const size_t T,
                                const float complex advances[restrict static T],
                                const size_t shift, const int down) {
    for (size_t it = 0; it < T; it++) {
        fwrite(&(int16_t) { lrintf(crealf(carrier) * 32767.0f) }, sizeof(int16_t), 1, stdout);
        carrier = renormalize(carrier * advances[((down ? T - it : it) + shift) % T]);
    }
    return carrier;
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

    char * bytes = NULL;
    size_t linecap = 0;
    const ssize_t B = getdelim(&bytes, &linecap, 0, stdin);
    if (B <= 0) exit(EXIT_FAILURE);

    const size_t D = (B * 8 + bits_per_sweep - 1) / bits_per_sweep;
    assert(D <= S);

    unsigned bits = 0;
    unsigned short bits_filled = 0;

    /* djb2 */
    unsigned hash = 5381;

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

    /* emit four unshifted upsweeps, with continuous carrier phase across sweeps */
    carrier = emit_sweep(carrier, T, advances, 0, 0);
    carrier = emit_sweep(carrier, T, advances, 0, 0);
    carrier = emit_sweep(carrier, T, advances, 0, 0);
    carrier = emit_sweep(carrier, T, advances, 0, 0);

    /* two unshifted downsweeps */
    carrier = emit_sweep(carrier, T, advances, 0, 1);
    carrier = emit_sweep(carrier, T, advances, 0, 1);

    /* one upsweep, circularly shifted to encode length of message in bytes */
    carrier = emit_sweep(carrier, T, advances, (B - 1) * L, 0);

    /* one shifted upsweep per data symbol */
    for (int ibyte = 0; ibyte < B || bits_filled; ) {
        while (bits_filled >= bits_per_sweep) {
            const unsigned symbol = bits & (S - 1U);
            carrier = emit_sweep(carrier, T, advances, symbol * L, 0);
            bits >>= bits_per_sweep;
            bits_filled -= bits_per_sweep;
        }
        if (B == ibyte && bits_filled) {
            /* emit an extra symbol to complete the last byte if there are leftover bits */
            const unsigned symbol = bits & (1U << bits_filled);
            carrier = emit_sweep(carrier, T, advances, symbol * L, 0);
            bits >>= bits_filled;
            bits_filled = 0;
        }
        while (bits_filled < 8 && ibyte < B) {
            const unsigned char byte = bytes[ibyte++];
            hash = hash * 33U ^ byte;
            bits |= byte << bits_filled;
            bits_filled += 8;
        }
    }

    /* one upsweep for checksum (lowest bits of djb2 of data symbols) */
    carrier = emit_sweep(carrier, T, advances, (hash & (S - 1U)) * L, 0);

    for (size_t ioffset = 0; ioffset < T; ioffset++)
        fwrite(&(int16_t) { 0 }, sizeof(int16_t), 1, stdout);

    free(advances);
}
