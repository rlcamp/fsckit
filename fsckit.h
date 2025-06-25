#include <stddef.h>
#include <stdint.h>
#include <complex.h>

/* standard posix definition of an iovec */
#ifndef _STRUCT_IOVEC
#define _STRUCT_IOVEC
struct iovec {
    char * iov_base;
    size_t iov_len;
};
#endif

float complex fsckit(void (* emit_sample_func)(void *, const int16_t), void * emit_sample_ctx,
                     const float amplitude, const float fs, const float fc, const float bw,
                     const unsigned bits_per_sweep, const unsigned interleave, float complex carrier,
                     const size_t V, const struct iovec iov[restrict static V]);
