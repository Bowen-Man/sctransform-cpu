"""Strict unbatched SCT v2 fitting primitives."""
from __future__ import annotations

import math
import numpy as np

from .model import (
    initial_overdispersion_r45,
    intercept_at_overdispersion_r45,
    matrixstats_refined_mean_r45,
    r45_offset_from_library,
)
from .r_glmgampoi import fit_dispersion_port_native
from .port_native import make_mu_native, reciprocal_native


def _local_dispersion_trend(expression: np.ndarray, alpha: np.ndarray) -> np.ndarray:
    """glmGamPoi-style Gaussian weighted running median for one fit batch."""
    n = len(alpha)
    if n == 0:
        raise ValueError("empty fit batch")
    npoints = min(n, max(100, int(round(0.1 * n))))
    half = npoints // 2
    order = np.argsort(expression, kind="stable")
    values = alpha[order]
    if 2 * half + 1 > n:
        weights = np.exp(-0.5 * np.linspace(-3, 3, n) ** 2)
        by_value = np.argsort(values, kind="stable")
        index = np.searchsorted(np.cumsum(weights[by_value]), weights.sum() / 2)
        return np.full(n, values[by_value[index]])
    weights = np.exp(-0.5 * np.linspace(-3, 3, 2 * half + 1) ** 2)
    result = np.empty(n)
    for index in range(half, n - half):
        window = values[index - half:index + half + 1]
        by_value = np.argsort(window, kind="stable")
        selected = np.searchsorted(np.cumsum(weights[by_value]), weights.sum() / 2)
        result[index] = window[by_value[selected]]
    result[:half] = result[half]
    result[n - half:] = result[n - half - 1]
    restored = np.empty(n)
    restored[order] = result
    return restored


def fit_batch(counts: np.ndarray, library: np.ndarray,
              all_cell_means: np.ndarray,
              all_cell_variances: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Fit one ordered, at-most-500-gene batch using the pinned R-compatible path."""
    counts = np.asarray(counts, dtype=np.float64)
    library = np.asarray(library, dtype=np.float64)
    if counts.ndim != 2 or counts.shape[1] != len(library):
        raise ValueError("counts must be gene x selected-cell")
    offset = r45_offset_from_library(library)
    initial_beta = np.empty(counts.shape[0])
    alpha = np.empty(counts.shape[0])
    expression = np.empty(counts.shape[0])

    for index, row in enumerate(counts):
        initial_alpha = initial_overdispersion_r45(row, offset)
        normalized = np.fromiter(
            (float(value) / math.exp(float(off)) for value, off in zip(row, offset)),
            dtype=np.float64,
            count=len(row),
        )
        beta0 = math.log(matrixstats_refined_mean_r45(normalized))
        beta = intercept_at_overdispersion_r45(row, offset, initial_alpha, beta0)
        initial_beta[index] = beta
        mu = make_mu_native(beta, offset)
        alpha[index] = fit_dispersion_port_native(row, mu, use_r45_mathlib=True)
        expression[index] = mu.mean()

    trend = _local_dispersion_trend(expression, alpha)
    beta = np.asarray([
        intercept_at_overdispersion_r45(row, offset, dispersion, start)
        for row, dispersion, start in zip(counts, trend, initial_beta)
    ])
    theta = reciprocal_native(alpha)
    naive_theta = all_cell_means ** 2 / (all_cell_variances - all_cell_means)
    theta[naive_theta / theta < 1e-3] = np.inf
    return beta, theta
