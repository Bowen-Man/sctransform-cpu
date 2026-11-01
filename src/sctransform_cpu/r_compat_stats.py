"""Experimental R 4.5 stats-compatible primitives for SCT fidelity work.

Algorithmic reference: R Core bandwidths.R, bandwidths.c and ksmooth.c,
R 4.5 branch, GPL-2-or-later; sctransform R/utils.R, GPL-3. See NOTICE.md and
docs/VALIDATION.md for attribution and the pinned validation limits.
Not intended as a replacement for R's general-purpose stats functions.
"""
from __future__ import annotations
import math
import numpy as np
from scipy import sparse

from .port_native import (
    bw_sj_native,
    ksmooth_normal_native,
    regularization_dispersion_native,
    regularization_theta_native,
    row_gmean_native,
)


def row_statistics_r45(counts):
    """Reference-order sparse mean, sample variance, and log geometric mean.

    For a canonical, sorted gene-by-cell sparse matrix, traverse each row's
    nonzero counts in cell order, corresponding to sctransform's CSC traversal.
    The final exponential uses the platform scalar ``exp`` matching the
    pinned Windows Rcpp build; ``expm1`` is not operation-order equivalent.
    """
    matrix = sparse.csr_matrix(counts, dtype=np.float64, copy=True)
    matrix.sum_duplicates(); matrix.eliminate_zeros(); matrix.sort_indices()
    rows, cols = matrix.shape
    if cols < 2: raise ValueError('at least two cells required')
    means=np.empty(rows); variances=np.empty(rows)
    for i in range(rows):
        values=matrix.data[matrix.indptr[i]:matrix.indptr[i+1]]
        total=0.
        for raw in values:
            value=float(raw)
            total += value
        mean=total/cols
        variance=0.
        for raw in values:
            variance += math.pow(float(raw)-mean,2)
        variance=(variance + math.pow(mean,2)*(cols-len(values)))/(cols-1)
        means[i]=mean
        variances[i]=variance
    _, log_gmeans = row_gmean_native(matrix)
    return means, variances, log_gmeans


def bw_sj_r45(values: np.ndarray, bins: int = 1000) -> float:
    """Pinned R 4.5.2 Sheather-Jones STE bandwidth."""
    return bw_sj_native(values, bins)


def outlier_r45(x: np.ndarray, y: np.ndarray, threshold: float = 10.) -> np.ndarray:
    """Double-shift bin MAD outlier mask for SCT v2's one-covariate case."""
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    if len(x) != len(y) or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("outlier input must be aligned and finite")
    width = np.ptp(x)*bw_sj_r45(x)/2
    if width <= 0: raise ValueError("outlier bin width is zero")
    epsilon = np.finfo(float).eps
    scores = []
    for shift in (0., .5):
        start = float(x.min()-10*epsilon-shift*width)
        stop = float(x.max()+width)
        breaks = start+np.arange(int(np.floor((stop-start)/width))+1)*width
        labels = np.searchsorted(breaks, x, side="left")-1
        if np.any(labels < 0) or np.any(labels >= len(breaks)-1):
            raise ValueError("outlier bin construction failed")
        score = np.empty(len(x))
        for group in np.unique(labels):
            mask = labels == group
            center = np.median(y[mask])
            mad = np.median(np.abs(y[mask]-center))*1.4826
            score[mask] = (y[mask]-center)/(mad+epsilon)
        scores.append(np.abs(score))
    return np.minimum(scores[0],scores[1]) > threshold


def ksmooth_normal_r45(x: np.ndarray, y: np.ndarray, query: np.ndarray, bandwidth: float) -> np.ndarray:
    """R ksmooth normal bandwidth scale and finite-support C loop."""
    return ksmooth_normal_native(x, y, query, bandwidth)


def regularize_v2_r45(log_gmean, fit_indices, intercept, theta, means, variances,
                       mean_filtered_library, bandwidth_adjust=3.0):
    """Candidate for the unbatched sctransform 0.4.3 v2 regularization path.

    No fallback approximation: unsupported/degenerate cases fail explicitly.
    Return the outlier mask in fit order to enable stage-by-stage comparisons.
    """
    log_gmean, means, variances = (np.asarray(v, dtype=np.float64)
                                  for v in (log_gmean, means, variances))
    fit_indices = np.asarray(fit_indices, dtype=int)
    intercept, theta = (np.asarray(v, dtype=np.float64) for v in (intercept, theta))
    x = log_gmean[fit_indices]
    dispersion = regularization_dispersion_native(x, theta)
    outliers = (outlier_r45(x, dispersion) | outlier_r45(x, intercept)
                | ~np.isfinite(theta))
    poisson = (variances <= means) | (means < 1e-3)
    keep = ~outliers & ~poisson[fit_indices]
    if keep.sum() < 2:
        raise ValueError('fewer than two NB fits remain for R-style regularization')
    bandwidth = bw_sj_r45(x[keep])*bandwidth_adjust
    query = np.clip(log_gmean, x[keep].min(), x[keep].max())
    b0 = ksmooth_normal_r45(x[keep], intercept[keep], query, bandwidth)
    od = ksmooth_normal_r45(x[keep], dispersion[keep], query, bandwidth)
    slope = ksmooth_normal_r45(
        x[keep], np.full(keep.sum(), math.log(10.0)), query, bandwidth)
    with np.errstate(divide='ignore', invalid='ignore'):
        fitted_theta = regularization_theta_native(log_gmean, od)
    fitted_theta[poisson] = np.inf
    b0[poisson] = np.log(means[poisson])-np.log(mean_filtered_library)
    slope[poisson] = math.log(10.0)
    return b0, fitted_theta, slope, bandwidth, outliers
