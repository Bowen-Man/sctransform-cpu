/* Narrow C++ arithmetic compatibility routines for the pinned Rcpp build.
 *
 * This deliberately mirrors the scalar expressions in sctransform's
 * src/utils.cpp.  It has no Python or R dependency at runtime; compiling as
 * C++ with the same Rtools toolchain is required for the strict Windows bit
 * contract.
 */
#include "math.h"

#ifdef _WIN32
# define SCT_EXPORT __declspec(dllexport)
#else
# define SCT_EXPORT __attribute__((visibility("default")))
#endif

extern "C" SCT_EXPORT int sct_get_fpu_cw()
{
    unsigned short cw;
    __asm__ __volatile__("fnstcw %0" : "=m" (cw));
    return cw;
}

extern "C" SCT_EXPORT void sct_use_r_fpu()
{
    unsigned short cw;
    __asm__ __volatile__("fnstcw %0" : "=m" (cw));
    cw = (unsigned short)((cw & ~0x0300u) | 0x0300u);
    __asm__ __volatile__("fldcw %0" : : "m" (cw));
}

extern "C" SCT_EXPORT void sct_row_gmean(
    const double data[], const int row_index[], int nnz, int rows, int cols,
    double gmean_out[], double log_gmean_out[])
{
    for (int row = 0; row < rows; ++row) {
        gmean_out[row] = 0.0;
    }
    for (int k = 0; k < nnz; ++k) {
        gmean_out[row_index[k]] += log(data[k] + 1.0);
    }
    for (int row = 0; row < rows; ++row) {
        double value = exp(gmean_out[row] / cols) - 1.0;
        gmean_out[row] = value;
        log_gmean_out[row] = log10(value);
    }
}

extern "C" SCT_EXPORT double sct_fit_beta_one_group(
    const double counts[], const double offsets[], int samples,
    double overdispersion, double beta)
{
    for (int iter = 0; iter < 100; ++iter) {
        double dl = 0.0;
        double ddl = 0.0;
        bool all_zero = true;
        for (int sample = 0; sample < samples; ++sample) {
            const double count = counts[sample];
            all_zero = all_zero && count == 0.0;
            const double mu = exp(beta + offsets[sample]);
            const double denom = 1.0 + mu * overdispersion;
            dl += (count - mu) / denom;
            ddl += mu * (1.0 + count * overdispersion) / denom / denom;
        }
        if (all_zero) {
            return -INFINITY;
        }
        const double step = dl / ddl;
        beta += step;
        if (fabs(step) < 1e-8 || isnan(beta)) {
            break;
        }
    }
    return beta;
}

extern "C" SCT_EXPORT void sct_make_mu(
    double beta, const double offsets[], int samples, double out[])
{
    for (int sample = 0; sample < samples; ++sample) {
        out[sample] = exp(beta + offsets[sample]);
    }
}

extern "C" SCT_EXPORT double sct_cpp_exp(double value)
{
    return exp(value);
}

