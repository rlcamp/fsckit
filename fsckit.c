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
        fwrite(&(float) { crealf(carrier) }, sizeof(float), 1, stdout);
        carrier = renormalize(carrier * advances[((down ? T - it : it) + shift) % T]);
    }
    return carrier;
}

int main(void) {
    /* sample rate */
    float fs = 8000;

    /* bandwidth of modulation, also sample rate if critically sampled */
    float bw = 250.0f;

    /* centre frequency */
    float fc = 2000.0f;

    /* this parameter also controls the spreading factor */
    unsigned bits_per_sweep = 5;

    /* number of unique measurable symbols is 2^bits_per_sweep */
    const size_t S = 1U << bits_per_sweep;

    const unsigned data_symbols[] = { 2, 4, 6, 0, 1 };
    const size_t D = sizeof(data_symbols) / sizeof(data_symbols[0]);
    assert(D <= S);

    /* djb2 */
    unsigned hash = 5381;
    for (size_t id = 0; id < D; id++)
        hash = hash * 33U ^ data_symbols[id];

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
        fwrite(&(float) { 0 }, sizeof(float), 1, stdout);

    /* emit four unshifted upsweeps, with continuous carrier phase across sweeps */
    carrier = emit_sweep(carrier, T, advances, 0, 0);
    carrier = emit_sweep(carrier, T, advances, 0, 0);
    carrier = emit_sweep(carrier, T, advances, 0, 0);
    carrier = emit_sweep(carrier, T, advances, 0, 0);

    /* two unshifted downsweeps */
    carrier = emit_sweep(carrier, T, advances, 0, 1);
    carrier = emit_sweep(carrier, T, advances, 0, 1);

    /* one upsweep, circularly shifted to encode length of message in symbols */
    carrier = emit_sweep(carrier, T, advances, (D - 1) * L, 0);

    /* one shifted upsweep per data symbol */
    for (size_t id = 0; id < D; id++)
        carrier = emit_sweep(carrier, T, advances, data_symbols[id] * L, 0);

    /* one upsweep for checksum (lowest bits of djb2 of data symbols) */
    carrier = emit_sweep(carrier, T, advances, (hash & (S - 1U)) * L, 0);

    for (size_t ioffset = 0; ioffset < T; ioffset++)
        fwrite(&(float) { 0 }, sizeof(float), 1, stdout);

    free(advances);
}
