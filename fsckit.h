#include <stddef.h>
#include <stdint.h>
#include <complex.h>

float complex fsckit(void (* emit_sample_func)(void *, const int16_t), void * emit_sample_ctx,
                     const float amplitude, const float fs, const float fc, const float bw,
                     const unsigned bits_per_sweep, const unsigned interleave, float complex carrier,
                     const size_t B, const unsigned char bytes[restrict static B]);
