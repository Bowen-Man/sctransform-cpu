/* Scalar glmGamPoi 1.22 dispersion callbacks for the strict Windows path.
 *
 * The formulas and operation order are derived from glmGamPoi's
 * src/overdispersion.cpp (GPL-3; portions adapted there from DESeq2 under
 * LGPL-3-or-later).  This narrow adapter supports the intercept-only design
 * used by sctransform v2.  The optimizer and data-flow implementation remain
 * in Python; these callbacks only prevent MSVC/Rtools arithmetic differences
 * from changing PORT control flow on nearly flat likelihoods.
 */
#include <math.h>
#include <stddef.h>

#ifdef _WIN32
# define SCT_EXPORT __declspec(dllexport)
#else
# define SCT_EXPORT __attribute__((visibility("default")))
#endif

extern double sct_r45_digamma(double x);
extern double sct_r45_trigamma(double x);

static const double cr_correction_factor = 0.99;

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
