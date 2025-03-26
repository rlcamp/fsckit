#include <stddef.h>
#include <stdint.h>

void fsckit(size_t (* get_bytes_func)(void *, const size_t Bmax, char buf[restrict static Bmax]),
            void * get_bytes_ctx,
            void (* emit_sample_func)(void *, const int16_t), void * emit_sample_ctx,
            const float amplitude, const float fs, const float fc, const float bw,
            const unsigned bits_per_sweep, const unsigned interleave);
