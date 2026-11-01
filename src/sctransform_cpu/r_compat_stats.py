"""Experimental R 4.5 stats-compatible primitives for SCT fidelity work.

Algorithmic reference: R Core bandwidths.R, bandwidths.c and ksmooth.c,
R 4.5 branch, GPL-2-or-later; sctransform R/utils.R, GPL-3. See
docs/STRICT_ALIGNMENT.md for attribution and the remaining validation limits.
Not intended as a replacement for R's general-purpose stats functions.
"""
from __future__ import annotations
import math
from decimal import Decimal, localcontext
import numpy as np
from scipy.optimize import brentq
from scipy import sparse


def row_statistics_r45(counts):
    """Reference-order sparse mean, sample variance, and log geometric mean.

    For a canonical, sorted gene-by-cell sparse matrix, traverse each row's
    nonzero counts in cell order, corresponding to sctransform's CSC traversal.
    Decimal exp is an experimental portability aid, not bitwise R guarantee.
    """
    matrix = sparse.csr_matrix(counts, dtype=np.float64, copy=True)
    matrix.sum_duplicates(); matrix.eliminate_zeros(); matrix.sort_indices()
    rows, cols = matrix.shape
    if cols < 2: raise ValueError('at least two cells required')
    means=np.empty(rows); variances=np.empty(rows); log_gmeans=np.empty(rows)
    with localcontext() as ctx:
        ctx.prec=50
        for i in range(rows):
            values=matrix.data[matrix.indptr[i]:matrix.indptr[i+1]]
            total=0.; log_total=0.
            for raw in values:
                value=float(raw)
                total += value
                log_total += math.log(value+1.)
            mean=total/cols
            variance=0.
            for raw in values:
                variance += math.pow(float(raw)-mean,2)
            variance=(variance + math.pow(mean,2)*(cols-len(values)))/(cols-1)
            exp_value=float(Decimal.from_float(log_total/cols).exp())
            gmean=exp_value-1.
            means[i]=mean
            variances[i]=variance
            log_gmeans[i]=math.log10(gmean) if gmean>0 else -math.inf
    return means, variances, log_gmeans


def bw_sj_r45(values: np.ndarray, bins: int = 1000) -> float:
    """R-style binned Sheather–Jones STE; numerical root may differ from R."""
    x = np.asarray(values, dtype=np.float64)
    n = len(x)
    if x.ndim != 1 or n < 2 or not np.all(np.isfinite(x)) or bins < 1:
        raise ValueError("bw.SJ requires at least two finite values")
    step = float(np.ptp(x)*1.01/bins)
    scale = min(float(np.std(x, ddof=1)), float(np.subtract(*np.percentile(x,[75,25]))/1.349))
    if step <= 0 or scale <= 0: raise ValueError("bw.SJ undefined for constant or zero-scale data")
    # R's pair counter bins absolute values toward zero (except its small-n
    # integer-overflow guard), then counts pairs by their bin-index distance.
    if n <= bins//2 and (x.min()/step < np.iinfo(np.int32).min or x.max()/step > np.iinfo(np.int32).max):
        indices = np.trunc((x-x.min())/step).astype(np.int64)
    else:
        indices = np.trunc(x/step).astype(np.int64)
        indices -= indices.min()
    if np.any(indices < 0) or np.any(indices >= bins): raise ValueError("bandwidth bin index out of range")
    histogram = np.bincount(indices, minlength=bins).astype(np.float64)
    pairs = np.correlate(histogram, histogram, mode="full")[bins-1:]
    pairs[0] = (pairs[0]-n)/2
    grid = np.arange(bins, dtype=np.float64)
    def functional(h, order):
        z2 = (grid*step/h)**2
        keep = z2 < 1000
        v = z2[keep]
        poly = v*v-6*v+3 if order == 4 else v*v*v-15*v*v+45*v-15
        diagonal = 3*n if order == 4 else -15*n
        return float((2*np.dot(pairs[keep], np.exp(-v/2)*poly)+diagonal)/
                     (n*(n-1)*np.sqrt(2*np.pi)*h**(order+1)))
    pilot4, pilot6 = 1.24*scale*n**(-1/7), 1.23*scale*n**(-1/9)
    fourth, sixth = functional(pilot4, 4), -functional(pilot6, 6)
    if not np.isfinite(sixth) or sixth <= 0 or not np.isfinite(fourth) or fourth <= 0:
        raise ValueError("bw.SJ pilot estimate invalid")
    c = 1/(2*np.sqrt(np.pi)*n)
    factor = 1.357*(fourth/sixth)**(1/7)
    upper = 1.144*scale*n**(-1/5)
    lower = .1*upper
    def equation(h):
        f = functional(factor*h**(5/7), 4)
        return (c/f)**(1/5)-h if f > 0 else np.nan
    for i in range(100):
        left, right = equation(lower), equation(upper)
        if np.isfinite(left) and np.isfinite(right) and left*right <= 0:
            return float(brentq(equation, lower, upper, xtol=.1*lower))
        if i % 2 == 0: upper *= 1.2
        else: lower /= 1.2
    raise ValueError("bw.SJ root not bracketed")


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
    """R ksmooth normal bandwidth scale and finite support, one output value."""
    x, y, query = (np.asarray(v, dtype=np.float64) for v in (x,y,query))
    if len(x) != len(y) or len(x)==0 or bandwidth<=0: raise ValueError("invalid smoother input")
    order = np.argsort(x, kind="stable"); xs, ys = x[order], y[order]
    sigma = bandwidth*0.3706506
    radius = 4*sigma
    result = np.empty(len(query), dtype=np.float64)
    for i, q in enumerate(query):
        left, right = np.searchsorted(xs, (q-radius, q+radius), side="left")
        # R includes a point exactly at the upper edge.
        right = np.searchsorted(xs, q+radius, side="right")
        w = np.exp(-.5*((xs[left:right]-q)/sigma)**2)
        result[i] = np.dot(w, ys[left:right])/w.sum() if len(w) else np.nan
    return result


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
    dispersion = np.log10(1 + np.power(10., x)/theta)
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
    with np.errstate(divide='ignore', invalid='ignore'):
        fitted_theta = np.power(10., log_gmean)/(np.power(10., od)-1)
    fitted_theta[poisson] = np.inf
    b0[poisson] = np.log(means[poisson])-np.log(mean_filtered_library)
    return b0, fitted_theta, bandwidth, outliers
