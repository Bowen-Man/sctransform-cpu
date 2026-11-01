/* Minimal standalone utilities for R 4.5.2's PORT Fortran sources.
 *
 * Derived from src/library/stats/src/port.c in R (GPL-2-or-later).
 * This file and binaries linked with portsrc.f are GPL-3.0-or-later.
 * The exported API is intentionally limited to one-dimensional nlminb
 * reverse communication used by the strict glmGamPoi compatibility path.
 */
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
# define SCT_EXPORT __declspec(dllexport)
#else
# define SCT_EXPORT __attribute__((visibility("default")))
#endif

extern void dv7dfl_(const int *alg, const int *lv, double v[]);
extern void drmngb_(double b[], double d[], double *fx, double g[],
                    int iv[], int *liv, int *lv, int *n, double v[],
                    double x[]);
extern void drmnhb_(double b[], double d[], double *fx, double g[],
                    double h[], int iv[], int *lh, int *liv, int *lv,
                    int *n, double v[], double x[]);

double d1mach_(const int *index)
{
    switch (*index) {
    case 1: return DBL_MIN;
    case 2: return DBL_MAX;
    case 3: return 0.5 * DBL_EPSILON;
    case 4: return DBL_EPSILON;
    case 5: return 0.301029995663981195213738894724493027;
    default: return 0.0;
    }
}

double dd7tpr_(int *n, const double x[], const double y[])
{
    double result = 0.0;
    for (int i = 0; i < *n; ++i) result += x[i] * y[i];
    return result;
}

void ditsum_(const double d[], const double g[], int iv[], const int *liv,
             const int *lv, const int *n, double v[], const double x[])
{
    (void)d; (void)g; (void)iv; (void)liv; (void)lv;
    (void)n; (void)v; (void)x;
}

void dn2cvp_(const int iv[], int *liv, int *lv, int *p, const double v[])
{ (void)iv; (void)liv; (void)lv; (void)p; (void)v; }
void dn2rdp_(const int iv[], int *liv, int *lv, int *n,
             const double rd[], const double v[])
{ (void)iv; (void)liv; (void)lv; (void)n; (void)rd; (void)v; }
void ds7cpr_(const double c[], const int iv[], int *l, int *liv)
{ (void)c; (void)iv; (void)l; (void)liv; }

void dv2axy_(int *n, double w[], const double *a,
             const double x[], const double y[])
{
    for (int i = 0; i < *n; ++i) w[i] = (*a) * x[i] + y[i];
}

double dv2nrm_(int *n, const double x[])
{
    double scale = 0.0, sumsq = 1.0;
    for (int i = 0; i < *n; ++i) {
        double value = fabs(x[i]);
        if (value == 0.0) continue;
        if (scale < value) {
            double ratio = scale / value;
            sumsq = 1.0 + sumsq * ratio * ratio;
            scale = value;
        } else {
            double ratio = value / scale;
            sumsq += ratio * ratio;
        }
    }
    return scale == 0.0 ? 0.0 : scale * sqrt(sumsq);
}

void dv7cpy_(int *n, double dest[], const double src[])
{ if (*n) memmove(dest, src, (size_t)*n * sizeof(double)); }
void dv7scl_(int *n, double dest[], const double *scale, const double src[])
{ for (int i = 0; i < *n; ++i) dest[i] = (*scale) * src[i]; }
void dv7scp_(int *n, double dest[], double *constant)
{ for (int i = 0; i < *n; ++i) dest[i] = *constant; }
void dv7swp_(int *n, double x[], double y[])
{
    for (int i = 0; i < *n; ++i) {
        double temporary = x[i]; x[i] = y[i]; y[i] = temporary;
    }
}

void dv7ipr_(int *n, const int permutation[], double x[])
{
    double *copy = (double *)malloc((size_t)*n * sizeof(double));
    for (int i = 0; i < *n; ++i) copy[i] = x[permutation[i] - 1];
    memcpy(x, copy, (size_t)*n * sizeof(double));
    free(copy);
}

void dv7prm_(int *n, const int permutation[], double x[])
{
    double *copy = (double *)malloc((size_t)*n * sizeof(double));
    for (int i = 0; i < *n; ++i) copy[permutation[i] - 1] = x[i];
    memcpy(x, copy, (size_t)*n * sizeof(double));
    free(copy);
}

void i7copy_(int *n, int dest[], const int src[])
{ if (*n) memmove(dest, src, (size_t)*n * sizeof(int)); }
void i7pnvr_(int *n, int x[], const int y[])
{ for (int i = 0; i < *n; ++i) x[y[i] - 1] = i + 1; }

void divset_(const int *alg, int iv[], const int *liv,
             const int *lv, double v[])
{
    static const int miniv[] = {0, 82, 59, 103, 103};
    static const int minv[] = {0, 98, 71, 101, 85};
    int algorithm = *alg;
    int alg1 = (algorithm - 1) % 2 + 1;
    if (*liv >= 21) iv[20] = 0;
    if (*liv >= 51) iv[50] = algorithm;
    if (algorithm < 1 || algorithm > 2) { iv[0] = 67; return; }
    if (*liv < miniv[algorithm]) { iv[0] = 15; return; }
    if (*lv < minv[algorithm]) { iv[0] = 16; return; }
    dv7dfl_(&alg1, lv, v);
    iv[0] = 12;
    iv[2] = 0;
    iv[43] = miniv[algorithm];
    iv[44] = minv[algorithm];
    iv[41] = minv[algorithm] + 1;
    iv[16] = 200;
    iv[17] = 150;
    iv[18] = 0;
    iv[19] = 1;
    iv[57] = miniv[algorithm] + 1;
    iv[21] = 0;
    iv[22] = 0;
    iv[3] = 0;
    iv[23] = 1;
    if (alg1 >= 2) {
        iv[15] = 0;
        iv[24] = 1;
        iv[51] = 0;
        iv[52] = 0;
        iv[49] = 25;
        iv[48] = 47;
        v[30] = 0.0;
    }
}

SCT_EXPORT int sct_port_state_size(void) { return 128; }

SCT_EXPORT void sct_port_initialize(int iv[], int liv, double v[], int lv,
                                    int max_iterations)
{
    int algorithm = 2;
    memset(iv, 0, (size_t)liv * sizeof(int));
    memset(v, 0, (size_t)lv * sizeof(double));
    divset_(&algorithm, iv, &liv, &lv, v);
    iv[17] = max_iterations;
}

SCT_EXPORT void sct_port_iterate_bfgs(double bounds[], double scale[],
                                      double *fx, double gradient[], int iv[],
                                      int liv, int lv, double v[], double x[])
{
    int n = 1;
    drmngb_(bounds, scale, fx, gradient, iv, &liv, &lv, &n, v, x);
}

SCT_EXPORT void sct_port_iterate_hessian(double bounds[], double scale[],
                                         double *fx, double gradient[],
                                         double hessian[], int iv[], int liv,
                                         int lv, double v[], double x[])
{
    int n = 1, lh = 1;
    drmnhb_(bounds, scale, fx, gradient, hessian, iv, &lh, &liv, &lv,
            &n, v, x);
}
