#include <stddef.h>
#include <stdint.h>

void unfsckit(int16_t * (* get_next_sample_func)(int16_t **, void *), void * get_ctx,
              void (* put_bytes_func)(const unsigned char *, const size_t, void *), void * put_ctx,
              const float sample_rate, const float f_carrier, const float bandwidth,
              const unsigned bits_per_sweep);
