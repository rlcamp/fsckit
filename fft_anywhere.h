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

 This file provides a header for the accompanying last-resort fft implementation, as well as any
 shims which use the same API to provide a common interface to additional FFT implementations. Any
 such implementations are covered by their original licenses.
 */

#include <complex.h>
#include <stddef.h>

/* Plan a forward complex-to-complex transform. T may be any seven-smooth number */
struct planned_forward_fft * plan_forward_fft_of_length(const size_t T);

/* Plan an inverse FFT of the given length. T may be any seven-smooth number */
struct planned_inverse_fft * plan_inverse_fft_of_length(const size_t T);

/* Plan a forward real-to-complex transform. T must be a seven-smooth multiple of 4 */
struct planned_real_fft * plan_real_fft_of_length(const size_t T);

/* Plan an inverse complex-to-real transform. T must be a seven-smooth multiple of 4 */
struct planned_real_inverse_fft * plan_real_inverse_fft_of_length(const size_t T);

/* Evaluate a single, foward, complex-to-complex transform. Input and output must not alias */
void fft_evaluate_forward(float complex * restrict const out, const float complex * restrict const in, const struct planned_forward_fft * const plan);

/* Evaluate a single, inverse, complex-to-complex transform. Input and output must not alias */
void fft_evaluate_inverse(float complex * restrict const out, const float complex * restrict const in, const struct planned_inverse_fft * const plan);

/* Evaluate a single, foward, real-to-complex transform. Out must be of length T/2, and is logically
 of length T / 2 + 1, with the real component of the Nyquist bin stored in the imaginary component
 of the DC bin. Input and output must not alias */
void fft_evaluate_real(float complex * restrict const out, const float * restrict const in, const struct planned_real_fft * const plan);

/* Evaluate a single, inverse, complex-to-real transform, with the DC and Nyquist bin stored as
 above. Input and output must not alias. Unlike the other three transform types, this also destroys
 its input */
void fft_evaluate_real_inverse(float * restrict const out, float complex * restrict const in, const struct planned_real_inverse_fft * const plan);

/* Destroy the above plans */
void destroy_planned_forward_fft(struct planned_forward_fft * plan);
void destroy_planned_inverse_fft(struct planned_inverse_fft * plan);
void destroy_planned_real_fft(struct planned_real_fft * plan);
void destroy_planned_real_inverse_fft(struct planned_real_inverse_fft * plan);

void fft_recursive_3(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan);
void fft_recursive_4(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan);
void fft_recursive_5(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan);
void fft_recursive_7(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan);
void fft_recursive_8(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan);

void fft_recursive_by_2(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan);
void fft_recursive_by_3(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan);
void fft_recursive_by_4(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan);
void fft_recursive_by_5(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan);
void fft_recursive_by_7(float complex * restrict const out, const float complex * restrict const in, const size_t istride, const struct planned_forward_fft * const plan);
