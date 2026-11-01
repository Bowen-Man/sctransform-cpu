"""R-order scalar primitives used by the strict SCT v2 fitting path."""
from __future__ import annotations

import math
import numpy as np


def r45_offset_from_library(library: np.ndarray) -> np.ndarray:
    """Reproduce ``log(10^log10_umi)`` used by sctransform 0.4.3.

    Although this is algebraically ``log(library)``, the intermediate base-10
    operations are observable at binary64 precision.  glmGamPoi subsequently
    evaluates means as ``exp(beta + offset)``, not ``exp(beta) * library``.
    """
    values = np.asarray(library, dtype=np.float64)
    if values.ndim != 1 or np.any(~np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("library sizes must be a positive finite vector")
    return np.fromiter((math.log(10.0 ** math.log10(float(value)))
                        for value in values), dtype=np.float64, count=len(values))


def matrixstats_refined_mean_r45(values) -> float:
    """Windows matrixStats ``rowMeans2(..., refine=TRUE)`` scalar order."""
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or len(values) == 0:
        raise ValueError("expected a nonempty vector")
    total = 0.0
    for value in values:
        total += float(value)
    average = total / len(values)
    residual = 0.0
    for value in values:
        residual += float(value) - average
    return average + residual / len(values)


def initial_overdispersion_r45(y: np.ndarray, offset: np.ndarray) -> float:
    """glmGamPoi 1.22 moment start for a one-row dense matrix.

    ``DelayedMatrixStats::rowVars`` uses a centered second pass here, while
    the exposure mean benefits from compensated accumulation.  Both details
    affect the last bits of the first intercept fit.
    """
    y = np.asarray(y, dtype=np.float64)
    offset = np.asarray(offset, dtype=np.float64)
    if y.ndim != 1 or offset.shape != y.shape or len(y) < 2:
        raise ValueError("expected aligned one-dimensional y and offset")
    mean = matrixstats_refined_mean_r45(y)
    centered = 0.0
    for value in y:
        delta = float(value) - mean
        centered += delta * delta
    variance = centered / (len(y) - 1)
    exposure_mean = math.fsum(math.exp(float(value)) for value in offset) / len(offset)
    if mean == 0.0:
        return 0.0
    inverse_exposure = 1.0 / exposure_mean
    result = (variance - inverse_exposure * mean) / (mean * mean)
    return 0.0 if not math.isfinite(result) or result < 0.0 else result


def intercept_at_overdispersion_r45(y: np.ndarray, offset: np.ndarray,
                                    alpha: float, initial: float) -> float:
    """Port glmGamPoi 1.22 ``fitBeta_one_group`` for one intercept.

    The scalar loop intentionally preserves the C++ accumulation and division
    order. ``alpha`` is glmGamPoi's overdispersion (``1/theta``).
    """
    y = np.asarray(y, dtype=np.float64)
    offset = np.asarray(offset, dtype=np.float64)
    if y.ndim != 1 or offset.shape != y.shape or len(y) == 0:
        raise ValueError("expected aligned one-dimensional y and offset")
    if not math.isfinite(alpha) or alpha < 0.0:
        raise ValueError("overdispersion must be finite and non-negative")
    beta = float(initial)
    for _ in range(100):
        dl = 0.0
        ddl = 0.0
        all_zero = True
        for count_raw, offset_raw in zip(y, offset):
            count = float(count_raw)
            all_zero = all_zero and count == 0.0
            mu = math.exp(beta + float(offset_raw))
            denominator = 1.0 + mu * alpha
            dl += (count - mu) / denominator
            ddl += mu * (1.0 + count * alpha) / denominator / denominator
        if all_zero:
            return -math.inf
        step = dl / ddl
        beta += step
        if abs(step) < 1e-8 or math.isnan(beta):
            break
    return beta
