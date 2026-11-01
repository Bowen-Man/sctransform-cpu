from pathlib import Path

import numpy as np
from scipy import sparse

from sctransform_cpu import SCTConfig, sctransform


def test_cpu_pipeline_is_deterministic_and_float64():
    rng = np.random.default_rng(11)
    library = rng.lognormal(4.5, 0.35, size=80)
    base = rng.lognormal(-4.0, 1.0, size=120)
    mean = base[:, None] * library[None, :]
    counts = rng.negative_binomial(2.0, 2.0 / (2.0 + mean)).astype(np.int32)
    config = SCTConfig(min_cells=2, output="full", features=40)
    first = sctransform(sparse.csr_matrix(counts), config=config)
    second = sctransform(sparse.csr_matrix(counts), config=config)
    assert first.residuals.dtype == np.float64
    assert np.array_equal(first.hvg_indices, second.hvg_indices)
    assert np.array_equal(first.residuals, second.residuals)
    assert np.array_equal(first.theta, second.theta)
