"""CPU dense gene-tile residuals and stable sample moments."""
from __future__ import annotations

import numpy as np


def cpu_tile(counts, library, intercept, theta, min_variance, clip, *, dtype=np.float32):
    mu = np.exp(intercept[:, None]) * library[None, :]
    variance = np.maximum(mu + mu * mu / theta[:, None], min_variance)
    residual = np.clip((counts - mu) / np.sqrt(variance), -clip, clip).astype(dtype)
    residual64 = residual.astype(np.float64)
    return residual, residual64.mean(axis=1), residual64.var(axis=1, ddof=1)
