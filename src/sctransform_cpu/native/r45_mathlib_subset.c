/* Narrow positive-argument R 4.5 Mathlib compatibility subset.
 *
 * Derived from R 4.5 src/nmath/{polygamma,gamma,lgamma,lgammacor,
 * chebyshev}.c, GPL-2-or-later.  It intentionally exports only the three
 * scalar functions needed by glmGamPoi's positive dispersion path.
 */
#include <float.h>
#include <math.h>
#include <stddef.h>

#if defined(_WIN32)
# define EXPORT __declspec(dllexport)
#else
# define EXPORT
#endif

static const double bvalues[] = {
  1.00000000000000000e+00, -5.00000000000000000e-01,
  1.66666666666666667e-01, -3.33333333333333333e-02,
  2.38095238095238095e-02, -3.33333333333333333e-02,
  7.57575757575757576e-02, -2.53113553113553114e-01,
  1.16666666666666667e+00, -7.09215686274509804e+00,
  5.49711779448621554e+01, -5.29124242424242424e+02,
  6.19212318840579710e+03, -8.65802531135531136e+04,
  1.42551716666666667e+06, -2.72982310678160920e+07,
  6.01580873900642368e+08, -1.51163157670921569e+10,
  4.29614643061166667e+11, -1.37116552050883328e+13,
  4.88332318973593167e+14, -1.92965793419400681e+16
};

static double dpsifn_positive(double x, int n)
{
  int i, k, nn, np, nx, fn;
  double arg, den, eps, fx, rxsq, s, t, ta, tk, tst, tt, t1, t2;
  double xdmln, xdmy, xinc, xm, xmin;
  double trm[23], trmr[101];
  double xln = log(x);
  double lrg = 1.0 / (2.0 * DBL_EPSILON);
  if (n == 0 && x * xln > lrg) return -xln;
  if (n >= 1 && x > n * lrg) return exp(-n * xln) / n;
  {
    const double r1m5 = 0.301029995663981195213738894724493027;
    const double r1m4 = DBL_EPSILON * 0.5;
    const double wdtol = r1m4 > 0.5e-18 ? r1m4 : 0.5e-18;
    const double elim = 2.302 * (1021 * r1m5 - 3.0);
    const double rln0 = r1m5 * 53.0;
    const double rln = rln0 < 18.06 ? rln0 : 18.06;
    double fln = (rln > 3.0 ? rln : 3.0) - 3.0;
    const double yint = 3.50 + 0.40 * fln;
    const double slope = 0.21 + fln * (0.0006038 * fln + 0.008677);
    nn = n;
    fn = nn;
    t = (fn + 1) * xln;
    if (fabs(t) > elim) return 0.0;
    if (x < wdtol) return pow(x, (double)-n-1.0);
    xm = yint + slope * fn;
    xmin = (int)xm + 1;
    if (n != 0) {
      xm = -2.302 * rln - fmin(0.0, xln);
      arg = fmin(0.0, xm / n);
      eps = exp(arg);
      xm = fabs(arg) < 1.0e-3 ? -arg : 1.0 - eps;
      fln = x * xm / eps;
      xm = xmin - x;
      if (xm > 7.0 && fln < 15.0) {
        int terms = (int)fln + 1;
        np = n + 1;
        t = exp(-(n + 1) * xln);
        s = t;
        den = x;
        for (i = 1; i <= terms; ++i) {
          den += 1.0;
          trm[i] = pow(den, (double)-np);
          s += trm[i];
        }
        return s;
      }
    }
    xdmy = x;
    xdmln = xln;
    xinc = 0.0;
    if (x < xmin) {
      nx = (int)x;
      xinc = xmin - nx;
      xdmy = x + xinc;
      xdmln = log(xdmy);
    }
    t = fn * xdmln;
    t1 = xdmln + xdmln;
    t2 = t + xdmln;
    if (fmax(fabs(t), fmax(fabs(t1), fabs(t2))) > elim) return 0.0;
    t = fn * xdmln;
    {
      const double wdtol2 = DBL_EPSILON * 0.5;
      double tss = exp(-t);
      tt = 0.5 / xdmy;
      t1 = tt;
      tst = wdtol2 * tt;
      if (nn != 0) t1 = tt + 1.0 / fn;
      rxsq = 1.0 / (xdmy * xdmy);
      ta = 0.5 * rxsq;
      t = (fn + 1) * ta;
      s = t * bvalues[2];
      if (fabs(s) >= tst) {
        tk = 2.0;
        for (k = 4; k <= 22; ++k) {
          t = t * ((tk + fn + 1)/(tk + 1.0))*((tk + fn)/(tk + 2.0))*rxsq;
          trm[k] = t * bvalues[k-1];
          if (fabs(trm[k]) < tst) break;
          s += trm[k];
          tk += 2.0;
        }
      }
      s = (s + t1) * tss;
      if (xinc != 0.0) {
        nx = (int)xinc;
        np = nn + 1;
        if (nn == 0) {
          for (i = 1; i <= nx; ++i) s += 1.0 / (x + (nx-i));
        } else {
          xm = xinc - 1.0;
          fx = x + xm;
          for (i = 1; i <= nx; ++i) {
            trmr[i] = pow(fx, (double)-np);
            s += trmr[i];
            xm -= 1.0;
            fx = x + xm;
          }
        }
      }
      return nn == 0 ? s - xdmln : s;
    }
  }
}

static double chebyshev_eval(double x, const double *a, int n)
{
  double b0 = 0.0, b1 = 0.0, b2 = 0.0, twox = x * 2.0;
  int i;
  for (i = 1; i <= n; ++i) {
    b2 = b1; b1 = b0; b0 = twox * b1 - b2 + a[n-i];
  }
  return (b0 - b2) * 0.5;
}

static const double gamcs[22] = {
 +.8571195590989331421920062399942e-2,+.4415381324841006757191315771652e-2,
 +.5685043681599363378632664588789e-1,-.4219835396418560501012500186624e-2,
 +.1326808181212460220584006796352e-2,-.1893024529798880432523947023886e-3,
 +.3606925327441245256578082217225e-4,-.6056761904460864218485548290365e-5,
 +.1055829546302283344731823509093e-5,-.1811967365542384048291855891166e-6,
 +.3117724964715322277790254593169e-7,-.5354219639019687140874081024347e-8,
 +.9193275519859588946887786825940e-9,-.1577941280288339761767423273953e-9,
 +.2707980622934954543266540433089e-10,-.4646818653825730144081661058933e-11,
 +.7973350192007419656460767175359e-12,-.1368078209830916025799499172309e-12,
 +.2347319486563800657233471771688e-13,-.4027432614949066932766570534699e-14,
 +.6910051747372100912138336975257e-15,-.1185584500221992907052387126192e-15
};

static const double algmcs[5] = {
 +.1666389480451863247205729650822e+0,-.1384948176067563840732986059135e-4,
 +.9810825646924729426157171547487e-8,-.1809129475572494194263306266719e-10,
 +.6221098041892605227126015543416e-13
};

static double gamma_positive(double x)
{
  int i, n = (int)x;
  double y = x - n, value;
  --n;
  value = chebyshev_eval(y * 2.0 - 1.0, gamcs, 22) + .9375;
  if (n == 0) return value;
  if (n < 0) {
    n = -n;
    for (i = 0; i < n; ++i) value /= (x + i);
    return value;
  }
  for (i = 1; i <= n; ++i) value *= (y + i);
  return value;
}

static double lgammacor_positive(double x)
{
  if (x < 94906265.62425156) {
    double tmp = 10.0 / x;
    return chebyshev_eval(tmp*tmp*2.0 - 1.0, algmcs, 5) / x;
  }
  return 1.0 / (x * 12.0);
}

EXPORT double sct_r45_digamma(double x)
{ return x > 0.0 ? -dpsifn_positive(x, 0) : NAN; }

EXPORT double sct_r45_trigamma(double x)
{ return x > 0.0 ? dpsifn_positive(x, 1) : NAN; }

/* Rtools45 uses mingw-w64 11.0.1.  Its x86 exp implementation is the x87
 * routine in mingw-w64-crt/math/x86/exp.def.h, licensed under ZPL-2.1.
 * The copy below is reduced to the finite binary64 range needed here and was
 * modified on 2026-11-01 to expose a project-specific symbol.  Keeping this
 * implementation local prevents a newer MinGW/UCRT exp from silently
 * changing the pinned R 4.5 numerical contract.
 */
static const long double r45_exp_c0 = 1.44268798828125L;
static const long double r45_exp_c1 = 7.05260771340735992468e-6L;

static long double r45_expl_internal(long double x)
{
  long double result = 0.0L;
  __asm__ volatile (
       "fldl2e\n\t"
       "fmul %%st(1),%%st\n\t"
       "subq $8, %%rsp\n"
       "fnstcw 4(%%rsp)\n"
       "movzwl 4(%%rsp), %%eax\n"
       "orb $12, %%ah\n"
       "movw %%ax, (%%rsp)\n"
       "fldcw (%%rsp)\n"
       "frndint\n\t"
       "fld %%st(1)\n\t"
       "frndint\n\t"
       "fldcw 4(%%rsp)\n"
       "addq $8, %%rsp\n"
       "fld %%st(1)\n\t"
       "fldt %2\n\t"
       "fld %%st(2)\n\t"
       "fmul %%st(1),%%st\n\t"
       "fsubp %%st,%%st(2)\n\t"
       "fld %%st(4)\n\t"
       "fsub %%st(3),%%st\n\t"
       "fmulp %%st,%%st(1)\n\t"
       "faddp %%st,%%st(1)\n\t"
       "fldt %3\n\t"
       "fmul %%st(4),%%st\n\t"
       "faddp %%st,%%st(1)\n\t"
       "f2xm1\n\t"
       "fld1\n\t"
       "faddp\n\t"
       "fstp %%st(1)\n\t"
       "fscale\n\t"
       "fstp %%st(1)\n\t"
       "fstp %%st(1)\n\t"
       : "=t" (result) : "0" (x), "m" (r45_exp_c0), "m" (r45_exp_c1)
       : "ax", "dx");
  return result;
}

EXPORT double sct_r45_exp(double x)
{
  if (isnan(x)) return x;
  if (x == INFINITY) return INFINITY;
  if (x == -INFINITY) return 0.0;
  if (x == 0.0) return 1.0;
  if (x > 709.782712893384) return INFINITY;
  if (x < -745.1332191019411) return 0.0;
  return (double)r45_expl_internal((long double)x);
}

EXPORT double sct_r45_lgamma(double x)
{
  const double ls2pi = 0.91893853320467274178032973640562;
  if (!(x > 0.0)) return NAN;
  if (x < 1e-306) return -log(x);
  if (x <= 10.0) return log(fabs(gamma_positive(x)));
  if (x > 1e17) return x * (log(x) - 1.0);
  if (x > 4934720.0) return ls2pi + (x - 0.5)*log(x) - x;
  return ls2pi + (x - 0.5)*log(x) - x + lgammacor_positive(x);
}
