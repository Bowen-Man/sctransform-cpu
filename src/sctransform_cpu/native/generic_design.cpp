/* General-design glmGamPoi primitives used by the covariate and batch path.
 *
 * The operation order follows glmGamPoi 1.22.0 beta_estimation.cpp,
 * fisher_scoring_steps.h, deviance.h and overdispersion.cpp.  Armadillo
 * RcppArmadillo 15.2.3-1's legacy Armadillo 14.6.3 headers plus the reference
 * BLAS/LAPACK backend preserve the pinned native linear-algebra semantics
 * without requiring R at runtime.
 */
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#ifdef _WIN32
# define SCT_EXPORT __declspec(dllexport)
#else
# define SCT_EXPORT __attribute__((visibility("default")))
#endif

static arma::mat copy_matrix(const double values[], int rows, int columns)
{
    arma::mat result(rows, columns);
    for (int row = 0; row < rows; ++row)
        for (int column = 0; column < columns; ++column)
            result(row, column) = values[(size_t)row * columns + column];
    return result;
}

extern "C" {
double dnrm2_(const int*, const double*, const int*);
void dscal_(const int*, const double*, double*, const int*);
double ddot_(const int*, const double*, const int*, const double*, const int*);
void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
void dgemm_(const char*, const char*, const int*, const int*, const int*,
            const double*, const double*, const int*, const double*, const int*,
            const double*, double*, const int*, arma::blas_len, arma::blas_len);
void dgesv_(const int*, const int*, double*, const int*, int*, double*,
            const int*, int*);
}

/* R's qr.default uses dqrdc2 rather than LAPACK.  This full-rank subset is
 * the exact loop from R 4.5's src/appl/dqrdc2.f, expressed against the same
 * reference BLAS.  Returning false on pivoting keeps the supported contract
 * explicit instead of silently relabelling a degenerate design. */
static bool r_qr(std::vector<double>& x, int n, int p,
                 std::vector<double>& qraux)
{
    const int one = 1;
    const double tolerance = 1e-7;
    std::vector<double> original(p), work(p);
    std::vector<int> pivot(p);
    for (int column = 0; column < p; ++column) {
        pivot[column] = column;
        qraux[column] = dnrm2_(&n, &x[(size_t)column * n], &one);
        original[column] = qraux[column];
        work[column] = qraux[column] == 0.0 ? 1.0 : qraux[column];
    }
    int rank_limit = p;
    int upper = std::min(n, p);
    for (int column = 0; column < upper; ++column) {
        while (column < rank_limit
               && qraux[column] < work[column] * tolerance) {
            for (int row = 0; row < n; ++row) {
                double saved = x[(size_t)column * n + row];
                for (int next = column + 1; next < p; ++next)
                    x[(size_t)(next - 1) * n + row] =
                        x[(size_t)next * n + row];
                x[(size_t)(p - 1) * n + row] = saved;
            }
            int saved_pivot = pivot[column];
            double saved_aux = qraux[column];
            double saved_original = original[column];
            double saved_work = work[column];
            for (int next = column + 1; next < p; ++next) {
                pivot[next - 1] = pivot[next];
                qraux[next - 1] = qraux[next];
                original[next - 1] = original[next];
                work[next - 1] = work[next];
            }
            pivot[p - 1] = saved_pivot;
            qraux[p - 1] = saved_aux;
            original[p - 1] = saved_original;
            work[p - 1] = saved_work;
            --rank_limit;
        }
        if (column == n - 1) continue;
        int length = n - column;
        double norm = dnrm2_(&length, &x[(size_t)column * n + column], &one);
        if (norm == 0.0) continue;
        if (x[(size_t)column * n + column] != 0.0)
            norm = std::copysign(norm, x[(size_t)column * n + column]);
        double inverse = 1.0 / norm;
        dscal_(&length, &inverse, &x[(size_t)column * n + column], &one);
        x[(size_t)column * n + column] += 1.0;
        for (int next = column + 1; next < p; ++next) {
            double value = -ddot_(&length,
                                  &x[(size_t)column * n + column], &one,
                                  &x[(size_t)next * n + column], &one)
                           / x[(size_t)column * n + column];
            daxpy_(&length, &value, &x[(size_t)column * n + column], &one,
                   &x[(size_t)next * n + column], &one);
            if (qraux[next] != 0.0) {
                double ratio = std::abs(x[(size_t)next * n + column])
                               / qraux[next];
                double update = std::max(1.0 - ratio * ratio, 0.0);
                if (std::abs(update) >= 1e-6) {
                    qraux[next] *= std::sqrt(update);
                } else {
                    int remaining = n - column - 1;
                    qraux[next] = dnrm2_(&remaining,
                                         &x[(size_t)next * n + column + 1],
                                         &one);
                    original[next] = qraux[next];
                }
            }
        }
        qraux[column] = x[(size_t)column * n + column];
        x[(size_t)column * n + column] = -norm;
    }
    int rank = std::min(rank_limit, n);
    if (rank != p) return false;
    for (int column = 0; column < p; ++column)
        if (pivot[column] != column) return false;
    return true;
}

static double gp_deviance(double y, double mu, double theta)
{
    if (theta < 1e-6) {
        if (y == 0.0) return 2.0 * mu;
        return std::max(2.0 * (y * std::log(y / mu) - (y - mu)), 0.0);
    }
    if (y == 0.0) return 2.0 / theta * std::log(1.0 + mu * theta);
    double s1 = y * std::log((mu + y * mu * theta) /
                             (y + y * mu * theta));
    double s2 = 1.0 / theta * std::log((1.0 + mu * theta) /
                                       (1.0 + y * theta));
    return std::max(-2.0 * (s1 - s2), 0.0);
}

static arma::vec calculate_mu(const arma::mat& design,
                              const arma::vec& beta,
                              const arma::vec& exp_offset)
{
    arma::vec mu = arma::exp(design * beta) % exp_offset;
    mu.transform([](double value) {
        return std::min(std::max(value, 1e-50), 1e50);
    });
    return mu;
}

extern "C" SCT_EXPORT int sct_general_initial_beta(
    const double counts[], int genes, int cells,
    const double design_values[], int coefficients,
    const double offsets[], int apply_ridge, double beta_out[])
{
    const int one = 1;
    const int qr_rows = cells + (apply_ridge ? coefficients : 0);
    const double ridge = 1e-10 / cells;
    std::vector<double> qr((size_t)qr_rows * coefficients, 0.0);
    for (int row = 0; row < cells; ++row)
        for (int column = 0; column < coefficients; ++column)
            qr[(size_t)column * qr_rows + row] =
                design_values[(size_t)row * coefficients + column];
    if (apply_ridge)
        for (int column = 0; column < coefficients; ++column)
            qr[(size_t)column * qr_rows + cells + column] = ridge;
    std::vector<double> qraux(coefficients);
    if (!r_qr(qr, qr_rows, coefficients, qraux)) return 1;

    std::vector<double> q((size_t)qr_rows * coefficients, 0.0);
    for (int column = 0; column < coefficients; ++column)
        q[(size_t)column * qr_rows + column] = 1.0;
    int reflectors = std::min(coefficients, qr_rows - 1);
    for (int column = 0; column < coefficients; ++column) {
        for (int reflector = reflectors - 1; reflector >= 0; --reflector) {
            if (qraux[reflector] == 0.0) continue;
            double saved = qr[(size_t)reflector * qr_rows + reflector];
            qr[(size_t)reflector * qr_rows + reflector] = qraux[reflector];
            int length = qr_rows - reflector;
            double value = -ddot_(&length,
                                  &qr[(size_t)reflector * qr_rows + reflector], &one,
                                  &q[(size_t)column * qr_rows + reflector], &one)
                           / qr[(size_t)reflector * qr_rows + reflector];
            daxpy_(&length, &value,
                   &qr[(size_t)reflector * qr_rows + reflector], &one,
                   &q[(size_t)column * qr_rows + reflector], &one);
            qr[(size_t)reflector * qr_rows + reflector] = saved;
        }
    }

    std::vector<double> responses((size_t)cells * genes);
    for (int gene = 0; gene < genes; ++gene) {
        for (int cell = 0; cell < cells; ++cell) {
            double normalized = counts[(size_t)gene * cells + cell]
                                / std::exp(offsets[cell]);
            responses[(size_t)gene * cells + cell] =
                std::log(normalized + 1.0);
        }
    }
    std::vector<double> solution((size_t)coefficients * genes);
    const char transpose = 'T', normal = 'N';
    const double one_double = 1.0, zero = 0.0;
    const arma::blas_len one_character = 1;
    dgemm_(&transpose, &normal, &coefficients, &genes, &cells,
           &one_double, q.data(), &qr_rows, responses.data(), &cells,
           &zero, solution.data(), &coefficients,
           one_character, one_character);

    std::vector<double> r((size_t)coefficients * coefficients, 0.0);
    for (int column = 0; column < coefficients; ++column)
        for (int row = 0; row <= column; ++row)
            r[(size_t)column * coefficients + row] =
                qr[(size_t)column * qr_rows + row];
    std::vector<int> pivots(coefficients);
    int info = 0;
    dgesv_(&coefficients, &genes, r.data(), &coefficients, pivots.data(),
           solution.data(), &coefficients, &info);
    if (info != 0) return 2;
    for (int gene = 0; gene < genes; ++gene)
        for (int coefficient = 0; coefficient < coefficients; ++coefficient)
            beta_out[(size_t)gene * coefficients + coefficient] =
                solution[(size_t)gene * coefficients + coefficient];
    return 0;
}

extern "C" SCT_EXPORT int sct_general_fisher(
    const double counts[], int genes, int cells,
    const double design_values[], int coefficients,
    const double offsets[], const double dispersions[],
    int apply_ridge, double beta_values[], double mu_out[])
{
    arma::mat design = copy_matrix(design_values, cells, coefficients);
    arma::vec exp_offset(cells);
    for (int cell = 0; cell < cells; ++cell)
        exp_offset(cell) = std::exp(offsets[cell]);

    for (int gene = 0; gene < genes; ++gene) {
        arma::vec counts_gene(cells), beta(coefficients);
        for (int cell = 0; cell < cells; ++cell)
            counts_gene(cell) = counts[(size_t)gene * cells + cell];
        for (int coefficient = 0; coefficient < coefficients; ++coefficient)
            beta(coefficient) = beta_values[(size_t)gene * coefficients + coefficient];

        double theta = dispersions[gene];
        const double ridge = 1e-10 / cells;
        arma::mat ridge_penalty = arma::eye(coefficients, coefficients) * ridge;
        arma::mat ridge_squared = ridge_penalty.t() * ridge_penalty;
        arma::vec ridge_target = arma::zeros(coefficients);
        arma::vec mu = calculate_mu(design, beta, exp_offset);
        double dev_old = 0.0;
        for (int cell = 0; cell < cells; ++cell)
            dev_old += gp_deviance(counts_gene(cell), mu(cell), theta);
        if (apply_ridge)
            dev_old += cells * arma::as_scalar(
                (beta - ridge_target).t() * ridge_squared
                * (beta - ridge_target));

        for (int iteration = 0; iteration < 1000; ++iteration) {
            arma::vec weight = mu / (1.0 + theta * mu);
            arma::vec weight_sqrt = arma::sqrt(weight);
            arma::mat q, r;
            arma::vec score;
            if (apply_ridge) {
                arma::mat ridge_helper = std::sqrt((double)cells) * ridge_penalty;
                arma::mat extended_design = arma::join_cols(design, ridge_helper);
                arma::vec extended_weight = arma::join_cols(
                    weight_sqrt, arma::ones(coefficients));
                arma::vec extended_residual = arma::join_cols(
                    (counts_gene - mu) / mu,
                    -ridge_helper * (beta - ridge_target));
                arma::mat weighted_design =
                    extended_design.each_col() % extended_weight;
                if (!arma::qr_econ(q, r, weighted_design)) return 2;
                score = (q.each_col() % extended_weight).t()
                        * extended_residual;
            } else {
                arma::mat weighted_design = design.each_col() % weight_sqrt;
                if (!arma::qr_econ(q, r, weighted_design)) return 2;
                score = (q.each_col() % weight_sqrt).t()
                        * ((counts_gene - mu) / mu);
            }
            arma::vec step = arma::solve(arma::trimatu(r), score);

            double speed = 1.0;
            int line_iteration = 0;
            beta += step;
            arma::vec old_mu = mu;
            double dev = 0.0;
            while (true) {
                mu = calculate_mu(design, beta, exp_offset);
                dev = 0.0;
                for (int cell = 0; cell < cells; ++cell)
                    dev += gp_deviance(counts_gene(cell), mu(cell), theta);
                if (apply_ridge)
                    dev += cells * arma::as_scalar(
                        (beta - ridge_target).t() * ridge_squared
                        * (beta - ridge_target));
                double convergence = std::abs(dev - dev_old) /
                                     (std::abs(dev) + 0.1);
                double relative_mu = arma::max(mu / old_mu);
                if ((dev < dev_old && relative_mu < 1e5)
                    || convergence < 1e-8) break;
                if (line_iteration >= 100) {
                    dev = std::numeric_limits<double>::quiet_NaN();
                    break;
                }
                speed /= 2.0;
                beta -= step * speed;
                ++line_iteration;
            }
            double convergence = std::abs(dev - dev_old) /
                                 (std::abs(dev) + 0.1);
            dev_old = dev;
            if (std::isnan(convergence)) return 3;
            if (convergence < 1e-8) break;
        }
        for (int coefficient = 0; coefficient < coefficients; ++coefficient)
            beta_values[(size_t)gene * coefficients + coefficient] =
                beta(coefficient);
        for (int cell = 0; cell < cells; ++cell)
            mu_out[(size_t)gene * cells + cell] = mu(cell);
    }
    return 0;
}

static void cr_matrices(const arma::mat& design, const arma::vec& mu,
                        double theta, arma::mat& b, arma::mat& db,
                        arma::mat& d2b)
{
    arma::vec weight = 1.0 / (1.0 / mu + theta);
    arma::vec derivative = -weight % weight;
    arma::vec second = -2.0 * derivative % weight;
    b = design.t() * (design.each_col() % weight);
    db = design.t() * (design.each_col() % derivative);
    d2b = design.t() * (design.each_col() % second);
}

static const double cr_factor = 0.99;

extern "C" SCT_EXPORT double sct_general_cr_logdet(
    const double design_values[], int cells, int coefficients,
    const double mu_values[], double theta)
{
    arma::mat design = copy_matrix(design_values, cells, coefficients);
    arma::vec mu(mu_values, cells);
    arma::mat b, db, d2b, l, u, p;
    cr_matrices(design, mu, theta, b, db, d2b);
    if (!arma::lu(l, u, p, b)) return NAN;
    double log_determinant = arma::sum(arma::log(l.diag()));
    arma::vec diagonal = u.diag();
    for (arma::uword index = 0; index < diagonal.n_elem; ++index) {
        double value = diagonal(index);
        log_determinant += std::log(value < 1e-50 ? 1e-50 : value);
    }
    return -0.5 * log_determinant * cr_factor;
}

extern "C" SCT_EXPORT void sct_general_cr_derivatives(
    const double design_values[], int cells, int coefficients,
    const double mu_values[], double theta, double out[])
{
    arma::mat design = copy_matrix(design_values, cells, coefficients);
    arma::vec mu(mu_values, cells);
    arma::mat b, db, d2b;
    cr_matrices(design, mu, theta, b, db, d2b);
    arma::mat inverse = arma::inv_sympd(
        b + arma::eye(coefficients, coefficients) * 1e-6);
    arma::mat product = inverse * db;
    double first = arma::trace(product);
    double second = first * first - arma::trace(product * product)
                    + arma::trace(inverse * d2b);
    out[0] = -0.5 * first * cr_factor;
    out[1] = (0.5 * first * first - 0.5 * second) * cr_factor;
}

extern "C" SCT_EXPORT double sct_general_cr_score_derivative(
    const double design_values[], int cells, int coefficients,
    const double mu_values[], double theta)
{
    arma::mat design = copy_matrix(design_values, cells, coefficients);
    arma::vec mu(mu_values, cells);
    arma::mat b, db, d2b;
    cr_matrices(design, mu, theta, b, db, d2b);
    arma::mat inverse = arma::inv_sympd(
        b + arma::eye(coefficients, coefficients) * 1e-6);
    return -0.5 * arma::trace(inverse * db) * cr_factor;
}

extern "C" SCT_EXPORT void sct_general_mu(
    const double beta_values[], int genes, int coefficients,
    const double design_values[], int cells, const double offsets[],
    double out[])
{
    arma::mat design = copy_matrix(design_values, cells, coefficients);
    arma::mat beta = copy_matrix(beta_values, genes, coefficients);
    arma::mat linear = beta * design.t();
    for (int gene = 0; gene < genes; ++gene)
        for (int cell = 0; cell < cells; ++cell) {
            double value = std::exp(linear(gene, cell) + offsets[cell]);
            out[(size_t)gene * cells + cell] =
                std::min(std::max(value, 1e-50), 1e50);
        }
}
