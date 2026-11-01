"""Unbatched SCT v2 corrected UMI candidate (official vst ordering)."""
from __future__ import annotations
import numpy as np


def corrected_counts_tile(counts, library, intercept, slope, theta, min_variance, *, scale_factor=None):
    """Reverse *unclipped* Pearson residuals at the median log-library size.

    sctransform::vst calls correct() before clipping rv$y; the target model
    variance in correct() does not use the residual minimum-variance floor.
    """
    counts = np.asarray(counts, dtype=np.float64)
    library = np.asarray(library, dtype=np.float64)
    if scale_factor is not None and (not np.isfinite(scale_factor) or scale_factor <= 0):
        raise ValueError('scale_factor must be finite and positive')
    log_umi = np.log10(library)
    target = np.median(log_umi) if scale_factor is None else np.log10(scale_factor)
    mu = np.exp(intercept[:, None] + slope[:, None]*log_umi[None, :])
    var = np.maximum(mu + mu*mu/theta[:, None], min_variance)
    raw_residual = (counts - mu)/np.sqrt(var)
    target_mu = np.exp(intercept[:, None] + slope[:, None]*target)
    target_var = target_mu + target_mu*target_mu/theta[:, None]
    corrected = np.rint(target_mu + raw_residual*np.sqrt(target_var))
    return np.maximum(corrected, 0.0)
