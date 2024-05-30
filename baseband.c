#include <complex.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

static float complex renormalize(const float complex x) {
    /* assuming x is already near unity, renormalize to unity w/o div or sqrt */
    const float magsquared = crealf(x) * crealf(x) + cimagf(x) * cimagf(x);
    return x * (3.0f - magsquared) * 0.5f;
}

static void butterworth_biquads(float num[][3], float den[][3], size_t P, float fs, float fc) {
    /* number of poles must be even */
    assert(!(P % 2));

    /* prewarp corner frequency for bilinear transform */
    const float wc = 2.0f * tanf((float)M_PI * fc / fs);

    /* each stage implements a conjugate pair of analog poles */
    for (size_t is = 0; is < P / 2; is++) {
        /* analog butterworth pole. the two poles for this stage are this and its conjugate */
        const float complex apole = wc * cexpf(I * (float)M_PI * (2.0f * (is + 1) + P - 1.0f) / (2 * P));

        /* each analog pole results in one digital pole and one digital zero at -1.0 */
        const float complex dpole = (2.0f - apole) / (2.0f + apole);

        /* polynomial coefficients for pair of digital zeros at -1 */
        num[is][0] = 1.0f;
        num[is][1] = 2.0f;
        num[is][2] = 1.0f;

        /* polynomial coefficients for conjugate pair of digital poles */
        den[is][0] = dpole * conjf(dpole);
        den[is][1] = -2.0f * crealf(dpole);
        den[is][2] = 1.0f;

        /* normalize the set of coefficients for unit gain */
        const float den_sum = den[is][0] + den[is][1] + den[is][2];
        const float den_scale = 1.0f / den[is][0], num_scale = den_scale * den_sum / 4.0f;
        for (size_t ik = 0; ik < 3; ik++) num[is][ik] *= num_scale;
        for (size_t ik = 0; ik < 3; ik++) den[is][ik] *= den_scale;
    }
}

static float complex cfilter(const float complex x, float complex vprev[2],
                             const float num[3], const float den[3]) {
    /* operate on complex input and output with real filter coefficients */
    const float complex v =          x - den[1] * vprev[0] - den[2] * vprev[1];
    const float complex y = num[0] * v + num[1] * vprev[0] + num[2] * vprev[1];

    vprev[1] = vprev[0];
    vprev[0] = v;

    return y;
}

int main(void) {
    const float sample_rate_out = 500;

    /* input arguments, all in cycles, samples, or symbols per second */
    const float sample_rate = 8000, f_carrier = 2000, bandwidth = 250;

    /* derived constants */
    float complex carrier = 1.0f;
    const float complex advance = cexpf(I * 2.0f * (float)M_PI * f_carrier / sample_rate);

    /* compute filter coefficients for eight-pole butterworth biquad cascade */
    float num[4][3], den[4][3];
    butterworth_biquads(num, den, 8, sample_rate, 0.5 * bandwidth);
    float complex vprev[4][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 } };

    const float input_samples_per_output_sample = sample_rate / sample_rate_out;
    float input_samples_since_output_sample = 0;

    float sample;
    while (fread(&sample, sizeof(float), 1, stdin)) {
        float complex filtered = sample * conjf(carrier);
        carrier = renormalize(carrier * advance);
        for (size_t is = 0; is < 4; is++)
            filtered = cfilter(filtered, vprev[is], num[is], den[is]);

        input_samples_since_output_sample++;
        if (input_samples_since_output_sample + 0.5f >= input_samples_per_output_sample) {
            input_samples_since_output_sample -= input_samples_per_output_sample;
            fwrite(&filtered, sizeof(float complex), 1, stdout);
        }
    }
}
