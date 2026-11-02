/* Scalar glmGamPoi 1.22 dispersion callbacks for strict platform paths.
 *
 * The formulas and operation order are derived from glmGamPoi's
 * src/overdispersion.cpp (GPL-3; portions adapted there from DESeq2 under
 * LGPL-3-or-later).  This narrow adapter supports the intercept-only design
 * used by sctransform v2. Python owns the high-level data flow; these routines
 * preserve the platform-specific arithmetic that controls PORT decisions and
 * the pinned residual/statistics outputs.
 */
#include <math.h>
#include <stddef.h>
#include <float.h>
#include <stdlib.h>

#ifdef _WIN32
# define SCT_EXPORT __declspec(dllexport)
#else
# define SCT_EXPORT __attribute__((visibility("default")))
#endif

extern double sct_r45_digamma(double x);
extern double sct_r45_trigamma(double x);
extern double sct_r45_lgamma(double x);

static const double cr_correction_factor = 0.99;

static double sct_platform_sqrt(double value)
{
#if defined(_WIN32) && defined(__GNUC__) && defined(__x86_64__)
    double result;
    __asm__ __volatile__("fldl %1; fsqrt; fstpl %0"
                         : "=m" (result) : "m" (value) : "st");
    return result;
#else
    return sqrt(value);
#endif
}

static double sct_platform_pow(double base, double exponent)
{
#ifdef _WIN32
    /* R's 64-bit Windows build defines USE_POWL_IN_R_POW. */
    return (double)powl((long double)base, (long double)exponent);
#else
    return pow(base, exponent);
#endif
}

static double sct_vector_lgamma(double value)
{
    return sct_r45_lgamma(value);
}

SCT_EXPORT void sct_log10_array(
    const double values[], int size, double out[])
{
    for (int i = 0; i < size; ++i) out[i] = log10(values[i]);
}

SCT_EXPORT void sct_reciprocal_array(
    const double values[], int size, double out[])
{
    for (int i = 0; i < size; ++i)
        out[i] = values[i] > 0.0 ? 1.0 / values[i] : INFINITY;
}

SCT_EXPORT void sct_regularization_dispersion(
    const double log_gmean[], const double theta[], int size, double out[])
{
    for (int i = 0; i < size; ++i) {
        double expression = sct_platform_pow(10.0, log_gmean[i]);
        out[i] = log10(1.0 + expression / theta[i]);
    }
}

SCT_EXPORT void sct_regularization_theta(
    const double log_gmean[], const double dispersion[], int size,
    double out[])
{
    for (int i = 0; i < size; ++i) {
        double expression = sct_platform_pow(10.0, log_gmean[i]);
        double denominator = sct_platform_pow(10.0, dispersion[i]) - 1.0;
        out[i] = expression / denominator;
    }
}

/* stats::var's scalar arithmetic is compiled as C in R itself.  Keeping this
 * adapter in a C translation unit (rather than the C++ compatibility unit)
 * is observable at the last bit for x87 long-double reductions. */
SCT_EXPORT void sct_count_mean_variance(
    const double values[], int size, double out[])
{
    long double total = 0.0L;
    for (int i = 0; i < size; ++i) total += values[i];
    long double extended_mean = total / size;
    double mean = (double)extended_mean;
    long double squared = 0.0L;
    for (int i = 0; i < size; ++i) {
        /* R stores the refined mean as binary64, then its x87 build keeps the
         * subtraction in an extended register before the explicit LDOUBLE
         * multiplication in cov.c. */
        long double difference = (long double)values[i] - (long double)mean;
        squared += difference * difference;
    }
    double variance = (double)(squared / (size - 1));
    out[0] = mean;
    out[1] = variance;
    volatile double numerator = variance - mean;
    volatile double denominator = mean * mean;
    volatile double start = numerator / denominator;
    out[2] = start;
    out[3] = log(start);
}

SCT_EXPORT void sct_ksmooth_normal(
    const double x[], const double y[], int n,
    const double query[], int query_size, double bandwidth, double out[])
{
    int imin = 0;
    double bw = bandwidth * 0.3706506;
    double cutoff = 4.0 * bw;
    while (imin < n && x[imin] < query[0] - cutoff) ++imin;
    for (int j = 0; j < query_size; ++j) {
        double num = 0.0, den = 0.0;
        double x0 = query[j];
        for (int i = imin; i < n; ++i) {
            if (x[i] < x0 - cutoff) {
                imin = i;
            } else {
                if (x[i] > x0 + cutoff) break;
                double distance = fabs(x[i] - x0) / bw;
                double weight = exp(-0.5 * distance * distance);
                num += weight * y[i];
                den += weight;
            }
        }
        out[j] = den > 0.0 ? num / den : NAN;
    }
}

typedef struct {
    int n;
    int bins;
    double d;
    double *counts;
    double c1;
    double alpha2;
} sct_bw_context;

static double sct_bw_phi4(const sct_bw_context *context, double h)
{
    double sum = 0.0;
    for (int i = 0; i < context->bins; ++i) {
        double delta = i * context->d / h;
        delta *= delta;
        if (delta >= 1000.0) break;
        double term = exp(-delta / 2.0)
                      * (delta * delta - 6.0 * delta + 3.0);
        sum += term * context->counts[i];
    }
    sum = 2.0 * sum + context->n * 3.0;
    return sum / ((double)context->n * (context->n - 1) * pow(h, 5.0))
           * 0.39894228040143267793994605993438;
}

static double sct_bw_phi6(const sct_bw_context *context, double h)
{
    double sum = 0.0;
    for (int i = 0; i < context->bins; ++i) {
        double delta = i * context->d / h;
        delta *= delta;
        if (delta >= 1000.0) break;
        double term = exp(-delta / 2.0)
                      * (delta * delta * delta - 15.0 * delta * delta
                         + 45.0 * delta - 15.0);
        sum += term * context->counts[i];
    }
    sum = 2.0 * sum - 15.0 * context->n;
    return sum / ((double)context->n * (context->n - 1) * pow(h, 7.0))
           * 0.39894228040143267793994605993438;
}

static double sct_bw_equation(double h, void *opaque)
{
    sct_bw_context *context = (sct_bw_context *)opaque;
    volatile double exponent = 5.0 / 7.0;
    volatile double pilot = context->alpha2 * pow(h, exponent);
    volatile double ratio = context->c1 / sct_bw_phi4(context, pilot);
    volatile double root = pow(ratio, 1.0 / 5.0);
    return root - h;
}

static double sct_zeroin2(double ax, double bx, double fa, double fb,
                          double (*function)(double, void *), void *info,
                          double tolerance, int maximum_iterations)
{
    double a = ax, b = bx, c = a, fc = fa;
    int remaining = maximum_iterations + 1;
    while (remaining--) {
        double previous_step = b - a;
        double p, q, new_step;
        if (fabs(fc) < fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }
        double actual_tolerance = 2.0 * DBL_EPSILON * fabs(b)
                                  + tolerance / 2.0;
        new_step = (c - b) / 2.0;
        if (fabs(new_step) <= actual_tolerance || fb == 0.0) return b;
        if (fabs(previous_step) >= actual_tolerance && fabs(fa) > fabs(fb)) {
            double t1, cb = c - b, t2;
            if (a == c) {
                t1 = fb / fa;
                p = cb * t1;
                q = 1.0 - t1;
            } else {
                q = fa / fc; t1 = fb / fc; t2 = fb / fa;
                p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
                q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
            }
            if (p > 0.0) q = -q; else p = -p;
            if (p < (0.75 * cb * q - fabs(actual_tolerance * q) / 2.0)
                    && p < fabs(previous_step * q / 2.0))
                new_step = p / q;
        }
        if (fabs(new_step) < actual_tolerance)
            new_step = new_step > 0.0 ? actual_tolerance : -actual_tolerance;
        a = b; fa = fb;
        b += new_step; fb = function(b, info);
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a; fc = fa;
        }
    }
    return b;
}

SCT_EXPORT double sct_bw_sj(const double sorted_x[], int n, int bins)
{
    double *counts = (double *)calloc((size_t)bins, sizeof(double));
    if (counts == NULL || n < 2 || bins < 1) {
        free(counts);
        return NAN;
    }
    double d = (sorted_x[n - 1] - sorted_x[0]) * 1.01 / bins;
    if (!(d > 0.0)) {
        free(counts);
        return NAN;
    }
    for (int i = 1; i < n; ++i) {
        int ii = (int)(sorted_x[i] / d);
        for (int j = 0; j < i; ++j) {
            int jj = (int)(sorted_x[j] / d);
            int distance = abs(ii - jj);
            if (distance < bins) counts[distance] += 1.0;
        }
    }

    double stats[4];
    sct_count_mean_variance(sorted_x, n, stats);
    double sd = sqrt(stats[1]);
    double q1_position = (n - 1) * 0.25;
    double q3_position = (n - 1) * 0.75;
    int q1_index = (int)floor(q1_position), q3_index = (int)floor(q3_position);
    double q1_fraction = q1_position - q1_index;
    double q3_fraction = q3_position - q3_index;
    double q1 = (1.0 - q1_fraction) * sorted_x[q1_index]
                + q1_fraction * sorted_x[q1_index + (q1_fraction > 0.0)];
    double q3 = (1.0 - q3_fraction) * sorted_x[q3_index]
                + q3_fraction * sorted_x[q3_index + (q3_fraction > 0.0)];
    double scale = fmin(sd, (q3 - q1) / 1.349);
    sct_bw_context context = {n, bins, d, counts, 0.0, 0.0};
    double a = 1.24 * scale * pow((double)n, -1.0 / 7.0);
    double b = 1.23 * scale * pow((double)n, -1.0 / 9.0);
    context.c1 = 1.0 / (2.0 * sqrt(3.141592653589793238462643383279502884)
                        * n);
    double td = -sct_bw_phi6(&context, b);
    context.alpha2 = 1.357 * pow(sct_bw_phi4(&context, a) / td, 1.0 / 7.0);
    double upper = 1.144 * scale * pow((double)n, -1.0 / 5.0);
    double lower = 0.1 * upper;
    double tolerance = 0.1 * lower;
    double left = sct_bw_equation(lower, &context);
    double right = sct_bw_equation(upper, &context);
    for (int attempt = 1; left * right > 0.0 && attempt <= 99; ++attempt) {
        if (attempt % 2) upper *= 1.2; else lower /= 1.2;
        left = sct_bw_equation(lower, &context);
        right = sct_bw_equation(upper, &context);
    }
    double result = sct_zeroin2(
        lower, upper, left, right, sct_bw_equation, &context, tolerance, 1000);
    free(counts);
    return result;
}

SCT_EXPORT void sct_pearson_residuals(
    const double counts[], int genes, int cells, const double log_umi[],
    const double intercept[], const double slope[], const double theta[],
    double minimum_variance, double clip, double out[],
    double mu_out[], double variance_out[], double denominator_out[])
{
    for (int gene = 0; gene < genes; ++gene) {
        for (int cell = 0; cell < cells; ++cell) {
            int index = gene * cells + cell;
            double linear = intercept[gene] + slope[gene] * log_umi[cell];
            double mu = exp(linear);
            double model_variance = mu + mu * mu / theta[gene];
            if (model_variance < minimum_variance)
                model_variance = minimum_variance;
            /* R materializes y-mu and sqrt(model_var) as separate vectors
             * before division.  Volatile stores preserve those binary64
             * boundaries instead of retaining an x87 extended numerator. */
            volatile double numerator = counts[index] - mu;
            volatile double denominator = sct_platform_sqrt(model_variance);
            double residual = numerator / denominator;
            if (residual < -clip) residual = -clip;
            if (residual > clip) residual = clip;
            out[index] = residual;
            if (mu_out != NULL) mu_out[index] = mu;
            if (variance_out != NULL) variance_out[index] = model_variance;
            if (denominator_out != NULL) denominator_out[index] = denominator;
        }
    }
}

SCT_EXPORT void sct_pearson_residuals_from_mu(
    const double counts[], const double mu[], int genes, int cells,
    const double theta[], double minimum_variance, double clip, double out[])
{
    for (int gene = 0; gene < genes; ++gene) {
        for (int cell = 0; cell < cells; ++cell) {
            int index = gene * cells + cell;
            double mean = mu[index];
            double model_variance = mean + mean * mean / theta[gene];
            if (model_variance < minimum_variance)
                model_variance = minimum_variance;
            volatile double numerator = counts[index] - mean;
            volatile double denominator = sct_platform_sqrt(model_variance);
            double residual = numerator / denominator;
            if (residual < -clip) residual = -clip;
            if (residual > clip) residual = clip;
            out[index] = residual;
        }
    }
}

/* Reproduce the two distinct reductions used by sctransform 0.4.3:
 *
 *   gene_attr$residual_mean     <- base::rowMeans(rv$y)
 *   gene_attr$residual_variance <- matrixStats::rowVars(rv$y)
 *
 * R's dense rowMeans accumulates in LDOUBLE (long double on the pinned
 * Windows/Rtools45 build).  matrixStats 1.5.0 rowVars instead uses binary64
 * throughout and performs its default refinement pass before the squared
 * deviations.  Keeping these paths separate is necessary for bit identity;
 * NumPy's pairwise/SIMD reductions deliberately use a different order. */
SCT_EXPORT void sct_residual_statistics(
    const double residual[], int genes, int cells,
    double mean_out[], double variance_out[])
{
    for (int gene = 0; gene < genes; ++gene) {
        const double *row = residual + (size_t)gene * (size_t)cells;

        long double mean_sum = 0.0L;
        for (int cell = 0; cell < cells; ++cell)
            mean_sum += (long double)row[cell];
        mean_out[gene] = (double)(mean_sum / (long double)cells);

        double sum = 0.0;
        for (int cell = 0; cell < cells; ++cell)
            sum += row[cell];
        double mu = sum / (double)cells;

        sum = 0.0;
        for (int cell = 0; cell < cells; ++cell)
            sum += row[cell] - mu;
        mu += sum / (double)cells;

        double sigma2 = 0.0;
        for (int cell = 0; cell < cells; ++cell) {
            double delta = row[cell] - mu;
            delta *= delta;
            sigma2 += delta;
        }
        variance_out[gene] = sigma2 / (double)(cells - 1);
    }
}

static void cox_reid_sums(const double mu[], int n, double theta,
                          double *b, double *db, double *d2b)
{
    double b_value = 0.0, db_value = 0.0, d2b_value = 0.0;
    for (int i = 0; i < n; ++i) {
        double weight = 1.0 / (1.0 / mu[i] + theta);
        double derivative = -1.0 * weight * weight;
        double second_derivative = -2.0 * derivative * weight;
        b_value += weight;
        db_value += derivative;
        d2b_value += second_derivative;
    }
    *b = b_value;
    *db = db_value;
    *d2b = d2b_value;
}

static void digamma_sums(const double y[], int n,
                         const double unique_counts[],
                         const double count_frequencies[], int table_size,
                         double theta_neg1,
                         double *digamma_term, double *trigamma_term,
                         double theta_neg2)
{
    double digamma_value = 0.0, trigamma_value = 0.0;
    if (table_size > 0) {
        for (int i = 0; i < table_size; ++i) {
            double argument = unique_counts[i] + theta_neg1;
            digamma_value += count_frequencies[i] * sct_r45_digamma(argument);
            trigamma_value += count_frequencies[i] * sct_r45_trigamma(argument);
        }
    } else {
        for (int i = 0; i < n; ++i) {
            double argument = y[i] + theta_neg1;
            digamma_value += sct_r45_digamma(argument);
            trigamma_value += sct_r45_trigamma(argument);
        }
    }
    trigamma_value *= theta_neg2;
    digamma_value -= n * sct_r45_digamma(theta_neg1);
    trigamma_value -= theta_neg2 * n * sct_r45_trigamma(theta_neg1);
    *digamma_term = digamma_value;
    *trigamma_term = trigamma_value;
}

SCT_EXPORT double sct_glmgampoi_score_derivative(
    const double y[], const double mu[], int n, double theta,
    int do_cox_reid_adjustment, const double unique_counts[],
    const double count_frequencies[], int table_size)
{
    double cr_term = 0.0, cr_term2 = 0.0;
    if (do_cox_reid_adjustment) {
        double b, db, d2b;
        cox_reid_sums(mu, n, theta, &b, &db, &d2b);
        double b_inv = 1.0 / (b + 1e-6);
        double d_i_db = b_inv * db;
        double ddetb = d_i_db;
        double d2detb = (ddetb * ddetb - d_i_db * d_i_db
                         + b_inv * d2b);
        cr_term = (0.5 * ddetb * ddetb - 0.5 * d2detb)
                  * cr_correction_factor;
        cr_term2 = -0.5 * ddetb * cr_correction_factor;
    }

    double theta_neg1 = 1.0 / theta;
    /* This matches the pinned Rtools build's R_pow_di(theta, -2). */
    double theta_neg2 = 1.0 / (theta * theta);
    double digamma_term, trigamma_term;
    digamma_sums(y, n, unique_counts, count_frequencies, table_size,
                 theta_neg1, &digamma_term, &trigamma_term, theta_neg2);

    double ll_part_1 = 0.0, ll_part_2 = 0.0;
    for (int i = 0; i < n; ++i) {
        ll_part_1 += log(1.0 + mu[i] * theta)
                     + (y[i] - mu[i]) / (mu[i] + theta_neg1);
        ll_part_2 += (mu[i] * mu[i] * theta + y[i])
                     / (1.0 + mu[i] * theta)
                     / (1.0 + mu[i] * theta);
    }
    double ll_part = -2.0 * theta_neg1 * (ll_part_1 - digamma_term)
                     + (ll_part_2 + trigamma_term);
    return ll_part + cr_term * (theta * theta)
           + (ll_part_1 - digamma_term) * theta_neg1
           + cr_term2 * theta;
}

SCT_EXPORT double sct_glmgampoi_score(
    const double y[], const double mu[], int n, double theta,
    int do_cox_reid_adjustment, const double unique_counts[],
    const double count_frequencies[], int table_size)
{
    double theta_neg1 = 1.0 / theta;
    double cr_term = 0.0;
    if (do_cox_reid_adjustment) {
        double b, db, ignored;
        cox_reid_sums(mu, n, theta, &b, &db, &ignored);
        double b_inv = 1.0 / (b + 1e-6);
        cr_term = -0.5 * (b_inv * db) * cr_correction_factor;
    }

    double digamma_term = 0.0;
    double max_y = 0.0, sum_y = 0.0, sum_prod_y = 0.0;
    if (table_size > 0) {
        for (int i = 0; i < table_size; ++i) {
            double value = unique_counts[i], frequency = count_frequencies[i];
            digamma_term += frequency * sct_r45_digamma(value + theta_neg1);
            sum_y += frequency * value;
            sum_prod_y += frequency * (value - 1.0) * value;
            max_y = fmax(max_y, value);
        }
    } else {
        for (int i = 0; i < n; ++i) {
            double value = y[i];
            digamma_term += sct_r45_digamma(value + theta_neg1);
            sum_y += value;
            sum_prod_y += (value - 1.0) * value;
            max_y = fmax(max_y, value);
        }
    }
    double correction = theta_neg1 > 1e5
        ? sum_prod_y / (2.0 * theta_neg1) : 0.0;
    if (max_y * 1e6 < theta_neg1) {
        digamma_term = sum_y - correction;
    } else {
        digamma_term -= n * sct_r45_digamma(theta_neg1);
        digamma_term *= theta_neg1;
        digamma_term = fmin(digamma_term, sum_y - correction);
    }

    double ll_part = 0.0;
    for (int i = 0; i < n; ++i) {
        double mu_theta = mu[i] * theta;
        if (mu_theta < 1e-10) {
            ll_part += mu_theta * mu_theta
                       * (1.0 / (1.0 + mu_theta) - 0.5);
        } else if (mu_theta < 1e-4) {
            double inverse = 1.0 / (1.0 + mu_theta);
            double upper = mu_theta * mu_theta * inverse;
            double lower = mu_theta * mu_theta * (inverse - 0.5);
            double suggested = log(1.0 + mu_theta)
                               - mu[i] / (mu[i] + theta_neg1);
            ll_part += fmax(fmin(suggested, upper), lower);
        } else {
            ll_part += log(1.0 + mu_theta)
                       - mu[i] / (mu[i] + theta_neg1);
        }
        ll_part += y[i] / (mu[i] + theta_neg1);
    }
    ll_part *= theta_neg1;
    return ll_part - digamma_term + cr_term * theta;
}

SCT_EXPORT double sct_glmgampoi_loglikelihood(
    const double y[], const double mu[], int n, double log_theta,
    int do_cox_reid_adjustment, const double unique_counts[],
    const double count_frequencies[], int table_size)
{
    double theta = exp(log_theta);
    double cr_term = 0.0;
    if (do_cox_reid_adjustment) {
        double b, ignored_db, ignored_d2b;
        cox_reid_sums(mu, n, theta, &b, &ignored_db, &ignored_d2b);
        cr_term = -0.5 * log(fmax(b, 1e-50)) * cr_correction_factor;
    }
    double theta_neg1 = 1.0 / theta;
    double lgamma_term = 0.0;
    if (table_size > 0) {
        for (int i = 0; i < table_size; ++i) {
            lgamma_term += count_frequencies[i]
                           * lgamma(unique_counts[i] + theta_neg1);
        }
    } else {
        for (int i = 0; i < n; ++i) {
            lgamma_term += sct_vector_lgamma(y[i] + theta_neg1);
        }
    }
    lgamma_term -= n * lgamma(theta_neg1);
    double ll_part = 0.0;
    for (int i = 0; i < n; ++i) {
        ll_part += (-y[i] - theta_neg1) * log(mu[i] + theta_neg1);
    }
    ll_part -= n * theta_neg1 * log(theta);
    return lgamma_term + ll_part + cr_term;
}

SCT_EXPORT double sct_glmgampoi_score_log(
    const double y[], const double mu[], int n, double log_theta,
    int do_cox_reid_adjustment, const double unique_counts[],
    const double count_frequencies[], int table_size)
{
    return sct_glmgampoi_score(
        y, mu, n, exp(log_theta), do_cox_reid_adjustment,
        unique_counts, count_frequencies, table_size);
}

SCT_EXPORT double sct_glmgampoi_score_derivative_log(
    const double y[], const double mu[], int n, double log_theta,
    int do_cox_reid_adjustment, const double unique_counts[],
    const double count_frequencies[], int table_size)
{
    return sct_glmgampoi_score_derivative(
        y, mu, n, exp(log_theta), do_cox_reid_adjustment,
        unique_counts, count_frequencies, table_size);
}

extern void sct_general_cr_derivatives(
    const double design_values[], int cells, int coefficients,
    const double mu_values[], double theta, double out[]);
extern double sct_general_cr_logdet(
    const double design_values[], int cells, int coefficients,
    const double mu_values[], double theta);
extern double sct_general_cr_score_derivative(
    const double design_values[], int cells, int coefficients,
    const double mu_values[], double theta);

SCT_EXPORT double sct_glmgampoi_general_loglikelihood(
    const double y[], const double mu[], int n, double log_theta,
    const double design[], int coefficients,
    const double unique_counts[], const double count_frequencies[],
    int table_size)
{
    double theta = exp(log_theta);
    double cr_term = sct_general_cr_logdet(
        design, n, coefficients, mu, theta);
    double theta_neg1 = 1.0 / theta;
    double lgamma_term = 0.0;
    if (table_size > 0) {
        for (int i = 0; i < table_size; ++i)
            lgamma_term += count_frequencies[i]
                           * lgamma(unique_counts[i] + theta_neg1);
    } else {
        for (int i = 0; i < n; ++i)
            lgamma_term += sct_vector_lgamma(y[i] + theta_neg1);
    }
    lgamma_term -= n * lgamma(theta_neg1);
    double ll_part = 0.0;
    for (int i = 0; i < n; ++i)
        ll_part += (-y[i] - theta_neg1) * log(mu[i] + theta_neg1);
    ll_part -= n * theta_neg1 * log(theta);
    return lgamma_term + ll_part + cr_term;
}

SCT_EXPORT double sct_glmgampoi_general_score(
    const double y[], const double mu[], int n, double log_theta,
    const double design[], int coefficients,
    const double unique_counts[], const double count_frequencies[],
    int table_size)
{
    double theta = exp(log_theta);
    double theta_neg1 = 1.0 / theta;
    double cr = sct_general_cr_score_derivative(
        design, n, coefficients, mu, theta);
    double digamma_term = 0.0;
    double max_y = 0.0, sum_y = 0.0, sum_prod_y = 0.0;
    if (table_size > 0) {
        for (int i = 0; i < table_size; ++i) {
            double value = unique_counts[i], frequency = count_frequencies[i];
            digamma_term += frequency * sct_r45_digamma(value + theta_neg1);
            sum_y += frequency * value;
            sum_prod_y += frequency * (value - 1.0) * value;
            max_y = fmax(max_y, value);
        }
    } else {
        for (int i = 0; i < n; ++i) {
            double value = y[i];
            digamma_term += sct_r45_digamma(value + theta_neg1);
            sum_y += value;
            sum_prod_y += (value - 1.0) * value;
            max_y = fmax(max_y, value);
        }
    }
    double correction = theta_neg1 > 1e5
        ? sum_prod_y / (2.0 * theta_neg1) : 0.0;
    if (max_y * 1e6 < theta_neg1) {
        digamma_term = sum_y - correction;
    } else {
        digamma_term -= n * sct_r45_digamma(theta_neg1);
        digamma_term *= theta_neg1;
        digamma_term = fmin(digamma_term, sum_y - correction);
    }
    double ll_part = 0.0;
    for (int i = 0; i < n; ++i) {
        double mu_theta = mu[i] * theta;
        if (mu_theta < 1e-10) {
            ll_part += mu_theta * mu_theta
                       * (1.0 / (1.0 + mu_theta) - 0.5);
        } else if (mu_theta < 1e-4) {
            double inverse = 1.0 / (1.0 + mu_theta);
            double upper = mu_theta * mu_theta * inverse;
            double lower = mu_theta * mu_theta * (inverse - 0.5);
            double suggested = log(1.0 + mu_theta)
                               - mu[i] / (mu[i] + theta_neg1);
            ll_part += fmax(fmin(suggested, upper), lower);
        } else {
            ll_part += log(1.0 + mu_theta)
                       - mu[i] / (mu[i] + theta_neg1);
        }
        ll_part += y[i] / (mu[i] + theta_neg1);
    }
    ll_part *= theta_neg1;
    return ll_part - digamma_term + cr * theta;
}

SCT_EXPORT double sct_glmgampoi_general_score_derivative_log(
    const double y[], const double mu[], int n, double log_theta,
    const double design[], int coefficients,
    const double unique_counts[], const double count_frequencies[],
    int table_size)
{
    double theta = exp(log_theta);
    double cr[2];
    sct_general_cr_derivatives(
        design, n, coefficients, mu, theta, cr);
    double theta_neg1 = 1.0 / theta;
    double theta_neg2 = 1.0 / (theta * theta);
    double digamma_term, trigamma_term;
    digamma_sums(y, n, unique_counts, count_frequencies, table_size,
                 theta_neg1, &digamma_term, &trigamma_term, theta_neg2);

    double ll_part_1 = 0.0, ll_part_2 = 0.0;
    for (int i = 0; i < n; ++i) {
        ll_part_1 += log(1.0 + mu[i] * theta)
                     + (y[i] - mu[i]) / (mu[i] + theta_neg1);
        ll_part_2 += (mu[i] * mu[i] * theta + y[i])
                     / (1.0 + mu[i] * theta)
                     / (1.0 + mu[i] * theta);
    }
    double ll_part = -2.0 * theta_neg1 * (ll_part_1 - digamma_term)
                     + (ll_part_2 + trigamma_term);
    return ll_part + cr[1] * (theta * theta)
           + (ll_part_1 - digamma_term) * theta_neg1
           + cr[0] * theta;
}
