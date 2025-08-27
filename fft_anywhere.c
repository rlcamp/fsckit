/*
 Copyright 2015-2025 Richard Campbell

 Permission to use, copy, modify, and/or distribute this software for any purpose with or without
 fee is hereby granted, provided that the above copyright notice and this permission notice appear
 in all copies.

 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
 SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
 AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
 THIS SOFTWARE.

 This is a good-enough FFT implementation which benchmarks within a respectable distance of the
 fastest known implementations on platforms of interest, and usually beats them on non-SIMD
 platforms, when compiled with a modern C compiler. This implementation should work where other FFT
 implementations do not, hence the name. If you are on a SIMD processor and need a permissive-
 licensed FFT implementation, jpommier/pffft is considerably faster, and you should prefer it to
 this or fftw.

 The core of the FFT is a set of functions for the DFTs/FFTs of size 3, 4, 5, 7, and 8. A second set
 of functions implements the Cooley-Tukey decomposition for T / S x S, for S = 2, 3, 4, 5, and 7.
 Each of these functions decomposes a length-T FFT into S FFTs of length T/S, followed by T/S FFTs
 of length S, where the latter are implemented using the relevant primitive DFT/FFT function. In
 this way, all FFTs of length T=2^a 3^b 5^c 7^d, for nonnegative integers a-d, may be computed.

 The twiddle factors are precomputed for a given FFT length and may be reused indefinitely and
 across threads. The transform is out-of-place, with the destination buffer used as scratch space.
 The inverse complex-to-real transform distorts its input, but the other three transforms do not.

 Heap allocation (using malloc) is performed only during fft planning, and is always freed
 in reverse order relative to allocation, allowing for simplified malloc implementations in
 an embedded firmware where malloc would otherwise not be needed.

 This code (and all high-performance C code using complex arithmetic) expects to be compiled with at
 least one of -ffinite-math-only, -fcx-limited-range, or -fcx-fortran-rules, in order to avoid a
 very significant slowdown due to the default semantics of complex multiplication. More modest
 speedups are obtained via (in descending order of benefit over risk ratio): "-fno-signed-zeros
 -fno-rounding-math -fexcess-precision=fast -fno-trapping-math -fno-math-errno -fassociative-math".
 In other words, this code is expected to perform fastest under -ffast-math semantics, but is almost
 as fast under "-ffast-math -fno-associative-math -fno-reciprocal-math", which should be palatable
 to a wider audience.
 */

/* needed for M_PI */
#define _XOPEN_SOURCE

#include "fft_anywhere.h"

#include <stdlib.h>
#include <math.h>

/* workaround for newlib and certain uncooperative combinations of compiler and libc */
#ifndef CMPLXF
#define CMPLXF __builtin_complex
#endif

struct planned_forward_fft {
    /* next function in recursive scheme */
    void (* nextfunc)(float complex * restrict, const float complex * restrict, size_t, const struct planned_forward_fft *);

    /* next c2c plan in the recursive scheme */
    struct planned_forward_fft * next;

    size_t T;

    /* twiddle factors, not used in the outermost c2c fft */
    float complex twiddles[];
};

static void dft3(float complex * restrict const out, const size_t stride, const float complex in0, const float complex in1, const float complex in2) {
    /* primitive for three-point discrete Fourier transform. this and the other primitives are inlined in several places */
    const float complex d1 = -0.5f - I * 0.866025404f;
    const float complex d2 = -0.5f + I * 0.866025404f;
    out[0 * stride] = in0 + in1 + in2;
    out[1 * stride] = in0 + d1 * in1 + d2 * in2;
    out[2 * stride] = in0 + d2 * in1 + d1 * in2;
}

static void dft5(float complex * restrict const out, const size_t stride, const float complex in0, const float complex in1, const float complex in2, const float complex in3, const float complex in4) {
    /* primitive for five-point discrete Fourier transform */
    const float complex d1 = +0.309016994f - I * 0.951056516f;
    const float complex d2 = -0.809016994f - I * 0.587785252f;
    const float complex d3 = -0.809016994f + I * 0.587785252f;
    const float complex d4 = +0.309016994f + I * 0.951056516f;
    out[0 * stride] = in0 +      in1 +      in2 +      in3 +      in4;
    out[1 * stride] = in0 + d1 * in1 + d2 * in2 + d3 * in3 + d4 * in4;
    out[2 * stride] = in0 + d2 * in1 + d4 * in2 + d1 * in3 + d3 * in4;
    out[3 * stride] = in0 + d3 * in1 + d1 * in2 + d4 * in3 + d2 * in4;
    out[4 * stride] = in0 + d4 * in1 + d3 * in2 + d2 * in3 + d1 * in4;
}

static void dft7(float complex * restrict const out, const size_t stride, const float complex in0, const float complex in1, const float complex in2, const float complex in3, const float complex in4, const float complex in5, const float complex in6) {
    /* primitive for seven-point discrete Fourier transform */
    const float complex d1 = +0.623489802f - I * 0.781831482f;
    const float complex d2 = -0.222520934f - I * 0.974927912f;
    const float complex d3 = -0.900968868f - I * 0.433883739f;
    const float complex d4 = -0.900968868f + I * 0.433883739f;
    const float complex d5 = -0.222520934f + I * 0.974927912f;
    const float complex d6 = +0.623489802f + I * 0.781831482f;
    out[0 * stride] = in0 +      in1 +      in2 +      in3 +      in4 +      in5 +      in6;
    out[1 * stride] = in0 + d1 * in1 + d2 * in2 + d3 * in3 + d4 * in4 + d5 * in5 + d6 * in6;
    out[2 * stride] = in0 + d2 * in1 + d4 * in2 + d6 * in3 + d1 * in4 + d3 * in5 + d5 * in6;
    out[3 * stride] = in0 + d3 * in1 + d6 * in2 + d2 * in3 + d5 * in4 + d1 * in5 + d4 * in6;
    out[4 * stride] = in0 + d4 * in1 + d1 * in2 + d5 * in3 + d2 * in4 + d6 * in5 + d3 * in6;
    out[5 * stride] = in0 + d5 * in1 + d3 * in2 + d1 * in3 + d6 * in4 + d4 * in5 + d2 * in6;
    out[6 * stride] = in0 + d6 * in1 + d5 * in2 + d4 * in3 + d3 * in4 + d2 * in5 + d1 * in6;
}

static void fft4(float complex * restrict const out, const size_t stride, const float complex in0, const float complex in1, const float complex in2, const float complex in3) {
    /* performs an fft of size 4 using four dft's of size 2, which results in 2/3 as many operations as a straight dft of size 4 */

    /* perform two dfts of size 2, one multiplied by a twiddle factor (a -90 degree phase shift) */
    const float complex scratch0 = in0 + in2;
    const float complex scratch1 = in0 - in2;
    const float complex scratch2 = in1 + in3;
    const float complex scratch3 = CMPLXF(cimagf(in1) - cimagf(in3), - crealf(in1) + crealf(in3));

    /* perform two more dfts of size 2 */
    out[0 * stride] = scratch0 + scratch2;
    out[1 * stride] = scratch1 + scratch3;
    out[2 * stride] = scratch0 - scratch2;
    out[3 * stride] = scratch1 - scratch3;
}

void fft_recursive_3(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan) { (void)plan;
    /* perform a three-point dft within the recursive framework */
    dft3(out, 1, in[0], in[istride], in[2 * istride]);
}

void fft_recursive_4(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan) { (void)plan;
    /* perform a four-point fft within the recursive framework */
    fft4(out, 1, in[0], in[istride], in[2 * istride], in[3 * istride]);
}

void fft_recursive_5(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan) { (void)plan;
    /* perform a five-point dft within the recursive framework */
    dft5(out, 1, in[0], in[istride], in[2 * istride], in[3 * istride], in[4 * istride]);
}

void fft_recursive_7(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan) { (void)plan;
    /* perform a five-point dft within the recursive framework */
    dft7(out, 1, in[0], in[istride], in[2 * istride], in[3 * istride], in[4 * istride], in[5 * istride], in[6 * istride]);
}

void fft_recursive_8(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan) { (void)plan;
    /* perform an eight-point fft within the recursive framework */
    const float complex in0 = in[0], in1 = in[istride], in2 = in[2 * istride], in3 = in[3 * istride], in4 = in[4 * istride], in5 = in[5 * istride], in6 = in[6 * istride], in7 = in[7 * istride];

    /* perform four dfts of size 2, two of which are multiplied by a twiddle factor (a -90 degree phase shift) */
    const float complex a0 = in0 + in4;
    const float complex a1 = in0 - in4;
    const float complex a2 = in2 + in6;
    const float complex a3 = CMPLXF(cimagf(in2) - cimagf(in6), crealf(in6) - crealf(in2));
    const float complex a4 = in1 + in5;
    const float complex a5 = in1 - in5;
    const float complex a6 = in3 + in7;
    const float complex a7 = CMPLXF(cimagf(in3) - cimagf(in7), crealf(in7) - crealf(in3));

    /* perform four more dfts of size 2 */
    const float complex c0 = a0 + a2;
    const float complex c1 = a1 + a3;
    const float complex c2 = a0 - a2;
    const float complex c3 = a1 - a3;
    const float complex c4 = a4 + a6;
    const float complex b5 = a5 + a7;
    const float complex b6 = a4 - a6;
    const float complex b7 = a5 - a7;

    /* apply final twiddle factors */
    const float complex c5 = CMPLXF(cimagf(b5) + crealf(b5),   cimagf(b5) - crealf(b5) ) * (float)M_SQRT1_2;
    const float complex c6 = CMPLXF(cimagf(b6), -crealf(b6));
    const float complex c7 = CMPLXF(cimagf(b7) - crealf(b7), -(crealf(b7) + cimagf(b7))) * (float)M_SQRT1_2;

    /* perform four dfts of length two */
    out[0] = c0 + c4;
    out[1] = c1 + c5;
    out[2] = c2 + c6;
    out[3] = c3 + c7;
    out[4] = c0 - c4;
    out[5] = c1 - c5;
    out[6] = c2 - c6;
    out[7] = c3 - c7;
}

void fft_recursive_by_2(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan) {
    const size_t T = plan->T;

    /* perform two ffts of length T / 2 */
    plan->nextfunc(out + 0 * (T / 2), in + 0 * istride, 2 * istride, plan->next);
    plan->nextfunc(out + 1 * (T / 2), in + 1 * istride, 2 * istride, plan->next);

    /* perform T / 2 dfts of length two, applying twiddle factors to all but the first */
    for (size_t it = 0; it < T / 2; it++) {
        const float complex tmp0 = out[it + 0 * (T / 2)];
        const float complex tmp1 = out[it + 1 * (T / 2)] * plan->twiddles[it];

        out[it + 0 * (T / 2)] = tmp0 + tmp1;
        out[it + 1 * (T / 2)] = tmp0 - tmp1;
    }
}

void fft_recursive_by_3(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft* const plan) {
    const size_t T = plan->T;

    /* perform three ffts of length T / 3 */
    plan->nextfunc(out + 0 * (T / 3), in + 0 * istride, 3 * istride, plan->next);
    plan->nextfunc(out + 1 * (T / 3), in + 1 * istride, 3 * istride, plan->next);
    plan->nextfunc(out + 2 * (T / 3), in + 2 * istride, 3 * istride, plan->next);

    /* perform T / 3 dfts of length three, applying twiddle factors to all but the first */
    for (size_t it = 0; it < T / 3; it++)
        dft3(out + it, T / 3,
             out[it + 0 * (T / 3)],
             out[it + 1 * (T / 3)] * plan->twiddles[2 * it + 0],
             out[it + 2 * (T / 3)] * plan->twiddles[2 * it + 1]);
}

void fft_recursive_by_4(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan) {
    const size_t T = plan->T;

    /* perform four ffts of length T / 4 */
    plan->nextfunc(out + 0 * (T / 4), in + 0 * istride, 4 * istride, plan->next);
    plan->nextfunc(out + 1 * (T / 4), in + 1 * istride, 4 * istride, plan->next);
    plan->nextfunc(out + 2 * (T / 4), in + 2 * istride, 4 * istride, plan->next);
    plan->nextfunc(out + 3 * (T / 4), in + 3 * istride, 4 * istride, plan->next);

    /* perform T / 4 ffts of length four, applying twiddle factors to all but the first */
    fft4(out + 0, T / 4, out[0], out[T / 4], out[2 * (T / 4)], out[3 * (T / 4)]);

    for (size_t it = 1; it < T / 4; it++)
        fft4(out + it, T / 4,
             out[it + 0 * (T / 4)],
             out[it + 1 * (T / 4)] * plan->twiddles[3 * it + 0],
             out[it + 2 * (T / 4)] * plan->twiddles[3 * it + 1],
             out[it + 3 * (T / 4)] * plan->twiddles[3 * it + 2]);
}

void fft_recursive_by_5(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan) {
    const size_t T = plan->T;

    /* perform five ffts of length T / 5 */
    plan->nextfunc(out + 0 * (T / 5), in + 0 * istride, 5 * istride, plan->next);
    plan->nextfunc(out + 1 * (T / 5), in + 1 * istride, 5 * istride, plan->next);
    plan->nextfunc(out + 2 * (T / 5), in + 2 * istride, 5 * istride, plan->next);
    plan->nextfunc(out + 3 * (T / 5), in + 3 * istride, 5 * istride, plan->next);
    plan->nextfunc(out + 4 * (T / 5), in + 4 * istride, 5 * istride, plan->next);

    /* perform T / 5 dfts of length five, applying twiddle factors to all but the first */
    for (size_t it = 0; it < T / 5; it++)
        dft5(out + it, T / 5,
             out[it + 0 * (T / 5)],
             out[it + 1 * (T / 5)] * plan->twiddles[4 * it + 0],
             out[it + 2 * (T / 5)] * plan->twiddles[4 * it + 1],
             out[it + 3 * (T / 5)] * plan->twiddles[4 * it + 2],
             out[it + 4 * (T / 5)] * plan->twiddles[4 * it + 3]);
}

void fft_recursive_by_7(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan) {
    const size_t T = plan->T;

    /* perform seven ffts of length T / 7 */
    plan->nextfunc(out + 0 * (T / 7), in + 0 * istride, 7 * istride, plan->next);
    plan->nextfunc(out + 1 * (T / 7), in + 1 * istride, 7 * istride, plan->next);
    plan->nextfunc(out + 2 * (T / 7), in + 2 * istride, 7 * istride, plan->next);
    plan->nextfunc(out + 3 * (T / 7), in + 3 * istride, 7 * istride, plan->next);
    plan->nextfunc(out + 4 * (T / 7), in + 4 * istride, 7 * istride, plan->next);
    plan->nextfunc(out + 5 * (T / 7), in + 5 * istride, 7 * istride, plan->next);
    plan->nextfunc(out + 6 * (T / 7), in + 6 * istride, 7 * istride, plan->next);

    /* perform T / 7 dfts of length seven, applying twiddle factors to all but the first */
    for (size_t it = 0; it < T / 7; it++)
        dft7(out + it, T / 7,
             out[it + 0 * (T / 7)],
             out[it + 1 * (T / 7)] * plan->twiddles[6 * it + 0],
             out[it + 2 * (T / 7)] * plan->twiddles[6 * it + 1],
             out[it + 3 * (T / 7)] * plan->twiddles[6 * it + 2],
             out[it + 4 * (T / 7)] * plan->twiddles[6 * it + 3],
             out[it + 5 * (T / 7)] * plan->twiddles[6 * it + 4],
             out[it + 6 * (T / 7)] * plan->twiddles[6 * it + 5]);
}

static float complex cosisinf(const float x) {
    return CMPLXF(cosf(x), sinf(x));
}

static void (* plan_recursive(struct planned_forward_fft ** plan_p, const size_t T))(float complex * restrict, const float complex * restrict, size_t, const struct planned_forward_fft *) {
    /* recursively plan a forward c2c fft for the given length, allocating and calculating all the
     necessary twiddle factors and branch conditions which will be encountered during execution of
     the fft, such that executing the fft only requires addition, multiplication, and following
     function pointers */

    if (T < 3) return NULL;

    if (3 == T) return fft_recursive_3;
    else if (4 == T) return fft_recursive_4;
    else if (5 == T) return fft_recursive_5;
    else if (7 == T) return fft_recursive_7;
    else if (8 == T) return fft_recursive_8;

    /* FFT size is not one of the primitive sizes, and must be divisible by a prime factor not larger than 7 */
    size_t S;

    if (6 == T || 10 == T || 14 == T) S = 2; /* special case so we don't end up at T = 2 */
    else if (!(T % 7)) S = 7;
    else if (!(T % 5)) S = 5;
    else if (!(T % 3)) S = 3;
    /* if T is not a power of 2, nothing we can do, calling code should check for this */
    else if (T & (T - 1)) return NULL;
    /* once reduced to powers of 2, try to get to repeatedly dividing by 4 and ending up at 8 */
    else S = T & (size_t)0xAAAAAAAAAAAAAAAA ? 4 : 2; /* S = 2 if T is a power of 4, else 4 */

    /* recursively plan the next fft size, if we can, before allocating the current one */
    struct planned_forward_fft * next = NULL;
    void (* nextfunc)(float complex * restrict, const float complex * restrict, size_t, const struct planned_forward_fft *) = plan_recursive(&next, T / S);
    if (!nextfunc) return NULL;

    /* only allocate once we are sure that we can plan the requested size */
    *plan_p = malloc(sizeof(struct planned_forward_fft) + sizeof(float complex) * (T / S) * (S - 1));
    **plan_p = (struct planned_forward_fft) { .T = T, .nextfunc = nextfunc, .next = next };

    for (size_t it = 0; it < T / S; it++)
        for (size_t is = 1; is < S; is++)
            (*plan_p)->twiddles[(is - 1) + (S - 1) * it] = cosisinf(-2.0f * it * is * (float)M_PI / T);

    return S == 7 ? fft_recursive_by_7 : S == 5 ? fft_recursive_by_5 : S == 4 ? fft_recursive_by_4 : S == 3 ? fft_recursive_by_3 : fft_recursive_by_2;
}

struct planned_forward_fft * plan_forward_fft_of_length(const size_t T) {
    struct planned_forward_fft * first = NULL;
    void (* nextfunc)(float complex * restrict, const float complex * restrict, size_t, const struct planned_forward_fft *) = plan_recursive(&first, T);
    if (!nextfunc) return NULL;

    /* only allocate once we are sure that we can plan the requested size */
    struct planned_forward_fft * plan = malloc(sizeof(*plan));
    *plan = (struct planned_forward_fft) { .T = T, .nextfunc = nextfunc, .next = first };
    return plan;
}

void destroy_planned_forward_fft(struct planned_forward_fft * plan) {
    if (!plan) return;
    struct planned_forward_fft * next = plan->next;
    free(plan);
    destroy_planned_forward_fft(next);
}

void destroy_planned_inverse_fft(struct planned_inverse_fft * plan) {
    destroy_planned_forward_fft((void *)plan);
}

void destroy_planned_real_fft(struct planned_real_fft * plan) {
    destroy_planned_forward_fft((void *)plan);
}

void destroy_planned_real_inverse_fft(struct planned_real_inverse_fft * plan) {
    destroy_planned_forward_fft((void *)plan);
}

void fft_evaluate_forward(float complex * restrict const out, const float complex * restrict const in, const struct planned_forward_fft * const plan) {
    plan->nextfunc(out, in, 1, plan->next);
}

struct planned_inverse_fft * plan_inverse_fft_of_length(const size_t T) {
    return (void *)plan_forward_fft_of_length(T);
}

void fft_evaluate_inverse(float complex * restrict const out, const float complex * restrict const in, const struct planned_inverse_fft * const iplan) {
    const struct planned_forward_fft * plan = (void *)iplan;
    const size_t T = plan->T;

    /* first compute the forward fft normally */
    plan->nextfunc(out, in, 1, plan->next);

    /* and then reverse the order of the outputs */
    for (size_t it = 1; it < T / 2; it++) {
        float complex tmp = out[it];
        out[it] = out[T - it];
        out[T - it] = tmp;
    }
}

struct planned_real_fft * plan_real_fft_of_length(const size_t T) {
    if ((T / 4U) * 4U != T) return NULL;

    /* recursively plan the next fft size, if we can, before allocating the current one */
    struct planned_forward_fft * next = NULL;
    void (* nextfunc)(float complex * restrict, const float complex * restrict, size_t, const struct planned_forward_fft *) = plan_recursive(&next, T / 2);
    if (!nextfunc) return NULL;

    struct planned_forward_fft * plan = malloc(sizeof(*plan) + sizeof(float complex) * (T / 4));
    *plan = (struct planned_forward_fft) { .nextfunc = nextfunc, .next = next, .T = T / 2 };

    for (size_t iw = 0; iw < T / 4; iw++)
        plan->twiddles[iw] = -I * cosisinf(-2.0f * (float)M_PI * iw / T);

    return (void *)plan;
}

struct planned_real_inverse_fft * plan_real_inverse_fft_of_length(const size_t T) {
    return (void *)plan_real_fft_of_length(T);
}

void fft_evaluate_real(float complex * restrict const out, const float * restrict const in, const struct planned_real_fft * const rplan) {
    const struct planned_forward_fft * const plan = (void *)rplan;
    const size_t Th = plan->T;

    plan->nextfunc(out, (void *)in, 1, plan->next);

    /* handle dc bin and nyquist bins. real component of nyquist bin is stored in imaginary component of dc bin */
    out[0] = CMPLXF(crealf(out[0]) + cimagf(out[0]), crealf(out[0]) - cimagf(out[0]));

    for (size_t iw = 1; iw < Th / 2; iw++) {
        /* this can probably be compactified more */
        const float complex a = out[iw], b = out[Th - iw];
        const float complex conj_b = conjf(b);
        const float complex tmpe = a + conj_b, tmpf = a - conj_b;
        const float complex tmpg = conjf(tmpe), tmph = -conjf(tmpf);
        const float complex tw = plan->twiddles[iw], conj_tw = conjf(tw);

        out[     iw] = 0.5f * (tmpe +      tw * tmpf);
        out[Th - iw] = 0.5f * (tmpg + conj_tw * tmph);
    }

    /* handle T/4 bin, for which the r2c twiddle factor is just -1.0 */
    out[Th / 2] = conjf(out[Th / 2]);
}

void fft_evaluate_real_inverse(float * restrict const out, float complex * restrict const in, const struct planned_real_inverse_fft * const iplan) {
    const struct planned_forward_fft * const plan = (void *)iplan;
    const size_t Th = plan->T, T = 2 * Th;

    in[Th / 2] = 2.0f * in[Th / 2];

    for (size_t iw = 1; iw < Th / 2; iw++) {
        const float complex a = in[iw], b = in[Th - iw];
        const float complex conj_b = conjf(b);
        const float complex tmpe = a + conj_b, tmpf = a - conj_b;
        const float complex tmpg = conjf(tmpe), tmph = conjf(tmpf);
        const float complex conj_tw = plan->twiddles[iw], tw = conjf(conj_tw);

        in[     iw] = tmpg + conj_tw * tmph;
        in[Th - iw] = tmpe -      tw * tmpf;
    }

    in[0] = CMPLXF(crealf(in[0]) + cimagf(in[0]), cimagf(in[0]) - crealf(in[0]));

    plan->nextfunc((void *)out, in, 1, plan->next);

    for (size_t it = 1; it < T; it += 2)
        out[it] = -out[it];
}
