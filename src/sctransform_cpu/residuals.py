"""CPU dense gene-tile residuals and stable sample moments."""
from __future__ import annotations

import numpy as np

from .port_native import (
    general_mu_native,
    pearson_residuals_from_mu_native,
    pearson_residuals_native,
    residual_statistics_native,
)


def cpu_tile(counts, library, intercept, slope, theta, min_variance, clip, *, dtype=np.float32):
    residual = pearson_residuals_native(
        counts, library, intercept, slope, theta, min_variance, clip
    ).astype(dtype, copy=False)
    mean, variance = residual_statistics_native(residual)
    return residual, mean, variance


def cpu_design_tile(counts, design, coefficients, theta, min_variance, clip):
    """Compute one dense gene tile for a general SCT design matrix."""
    mu = general_mu_native(coefficients, design, np.zeros(len(design)))
    residual = pearson_residuals_from_mu_native(
        counts, mu, theta, min_variance, clip)
    mean, variance = residual_statistics_native(residual)
    return residual, mean, variance
