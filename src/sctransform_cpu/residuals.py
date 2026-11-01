"""CPU dense gene-tile residuals and stable sample moments."""
from __future__ import annotations

import numpy as np

from .port_native import pearson_residuals_native, residual_statistics_native


def cpu_tile(counts, library, intercept, slope, theta, min_variance, clip, *, dtype=np.float32):
    residual = pearson_residuals_native(
        counts, library, intercept, slope, theta, min_variance, clip
    ).astype(dtype, copy=False)
    mean, variance = residual_statistics_native(residual)
    return residual, mean, variance
