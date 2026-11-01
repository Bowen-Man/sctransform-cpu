"""Minimal CPU-only implementation of the validated unbatched SCT v2 path."""
from __future__ import annotations

from dataclasses import dataclass, field
import time

import numpy as np
from scipy import sparse

from .corrected_counts import corrected_counts_tile
from .fit import fit_batch
from .r_compat_stats import regularize_v2_r45, row_statistics_r45
from .residuals import cpu_tile


@dataclass(frozen=True)
class SCTConfig:
    min_cells: int = 5
    fit_cells: int | None = None
    fit_genes: int | None = None
    features: int = 3000
    output: str = "hvg"
    tile_genes: int = 256
    bandwidth_adjust: float = 3.0
    return_corrected_umi: bool = False
    corrected_scale_factor: float | None = None


@dataclass
class SCTResult:
    genes: np.ndarray
    cells: np.ndarray
    intercept: np.ndarray
    theta: np.ndarray
    slope: np.ndarray
    residual_mean: np.ndarray
    residual_variance: np.ndarray
    hvg_indices: np.ndarray
    residual_gene_indices: np.ndarray
    residuals: np.ndarray
    fit_gene_indices: np.ndarray
    fit_cell_indices: np.ndarray
    fit_intercept: np.ndarray
    fit_theta: np.ndarray
    min_variance: float
    bandwidth: float
    timing: dict[str, float] = field(default_factory=dict)
    diagnostics: dict[str, object] = field(default_factory=dict)
    corrected_umi: sparse.csr_matrix | None = None
    semantics: str = (
        "platform-pinned Windows x64 or Linux x86_64 R 4.5.2 / "
        "sctransform 0.4.3 / "
        "glmGamPoi 1.22.0 / matrixStats 1.5.0; direct unbatched v2 path only"
    )


def _names(values, size: int, prefix: str) -> np.ndarray:
    result = np.asarray(
        values if values is not None else [f"{prefix}-{index}" for index in range(size)],
        dtype=str,
    )
    if len(result) != size or len(set(result)) != size:
        raise ValueError(f"{prefix} names must match the matrix and be unique")
    return result


def sctransform(counts, *, gene_names=None, cell_names=None,
                config: SCTConfig | None = None,
                fit_cell_indices=None, fit_gene_names=None) -> SCTResult:
    """Transform raw UMI counts in gene x cell orientation.

    Exact comparison with R requires the same ordered fit cells and genes.
    This release deliberately refuses implicit subsampling; pass a shared R
    selection through ``fit_cell_indices`` and ``fit_gene_names``.
    """
    cfg = config or SCTConfig()
    if cfg.min_cells < 1 or cfg.features < 1 or cfg.tile_genes < 1:
        raise ValueError("min_cells, features and tile_genes must be positive")
    if cfg.bandwidth_adjust <= 0 or cfg.output not in {"hvg", "full"}:
        raise ValueError("invalid bandwidth_adjust or output")
    if any(value is not None and value < 1 for value in (cfg.fit_cells, cfg.fit_genes)):
        raise ValueError("fit_cells and fit_genes must be positive or None")

    started = time.perf_counter()
    matrix = sparse.csr_matrix(counts, dtype=np.float64, copy=True)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    matrix.sort_indices()
    if matrix.ndim != 2 or matrix.shape[0] < 1 or matrix.shape[1] < 2:
        raise ValueError("expected a nonempty gene x cell matrix with at least two cells")
    if (np.any(~np.isfinite(matrix.data)) or np.any(matrix.data < 0)
            or np.any(matrix.data != np.floor(matrix.data))):
        raise ValueError("input must contain finite non-negative integer raw UMI counts")
    genes = _names(gene_names, matrix.shape[0], "gene")
    cells = _names(cell_names, matrix.shape[1], "cell")
    library = np.asarray(matrix.sum(axis=0)).ravel()
    keep = np.asarray(matrix.getnnz(axis=1)).ravel() >= cfg.min_cells
    matrix, genes = matrix[keep].tocsr(), genes[keep]
    if matrix.shape[0] < 2 or np.any(library <= 0):
        raise ValueError("too few retained genes or a zero-library cell")
    retained_mean_library = float(np.asarray(matrix.sum(axis=0)).mean())
    means, variances, log_gmean = row_statistics_r45(matrix)
    timing = {"input_statistics": time.perf_counter() - started}

    started = time.perf_counter()
    if fit_cell_indices is None:
        requested_cells = matrix.shape[1] if cfg.fit_cells is None else cfg.fit_cells
        if requested_cells < matrix.shape[1]:
            raise ValueError("implicit cell subsampling is disabled; provide fit_cell_indices")
        selected_cells = np.arange(matrix.shape[1])
    else:
        selected_cells = np.asarray(fit_cell_indices)
        if (selected_cells.ndim != 1
                or not np.issubdtype(selected_cells.dtype, np.integer)
                or len(selected_cells) == 0
                or len(np.unique(selected_cells)) != len(selected_cells)
                or np.any(selected_cells < 0)
                or np.any(selected_cells >= matrix.shape[1])):
            raise ValueError("fit_cell_indices must contain distinct valid integer indices")
    sampled = matrix[:, selected_cells]
    eligible = np.flatnonzero(
        (np.asarray(sampled.getnnz(axis=1)).ravel() >= cfg.min_cells)
        & (variances > means)
    )
    if fit_gene_names is not None:
        requested = np.asarray(fit_gene_names, dtype=str)
        lookup = {name: index for index, name in enumerate(genes)}
        if (requested.ndim != 1 or len(set(requested)) != len(requested)
                or any(name not in lookup for name in requested)):
            raise ValueError("fit_gene_names must contain distinct retained genes")
        chosen = np.asarray([lookup[name] for name in requested], dtype=int)
        if not np.isin(chosen, eligible).all():
            raise ValueError("fit_gene_names must satisfy v2 fitting eligibility")
        eligible = chosen
    elif cfg.fit_genes is not None and cfg.fit_genes < len(eligible):
        raise ValueError("implicit gene subsampling is disabled; provide fit_gene_names")
    if len(eligible) < 2:
        raise ValueError("fewer than two overdispersed genes are available for fitting")
    timing["selection"] = time.perf_counter() - started

    started = time.perf_counter()
    fit_intercept = np.empty(len(eligible))
    fit_theta = np.empty(len(eligible))
    for offset in range(0, len(eligible), 500):
        indices = eligible[offset:offset + 500]
        fit_intercept[offset:offset + len(indices)], fit_theta[offset:offset + len(indices)] = fit_batch(
            sampled[indices].toarray(),
            library[selected_cells],
            means[indices],
            variances[indices],
        )
    timing["fit"] = time.perf_counter() - started

    started = time.perf_counter()
    intercept, theta, slope, bandwidth, outliers = regularize_v2_r45(
        log_gmean,
        eligible,
        fit_intercept,
        fit_theta,
        means,
        variances,
        retained_mean_library,
        cfg.bandwidth_adjust,
    )
    timing["regularization"] = time.perf_counter() - started

    started = time.perf_counter()
    clip = float(np.sqrt(matrix.shape[1]))
    min_variance = float((np.median(matrix.data) / 5.0) ** 2)
    residual_mean = np.empty(len(genes))
    residual_variance = np.empty(len(genes))
    full = np.empty(matrix.shape, dtype=np.float64) if cfg.output == "full" else None
    corrected = [] if cfg.return_corrected_umi else None
    for offset in range(0, len(genes), cfg.tile_genes):
        indices = np.arange(offset, min(offset + cfg.tile_genes, len(genes)))
        dense = matrix[indices].toarray()
        residual, mean, variance = cpu_tile(
            dense, library, intercept[indices], slope[indices], theta[indices], min_variance, clip,
            dtype=np.float64,
        )
        residual_mean[indices], residual_variance[indices] = mean, variance
        if full is not None:
            full[indices] = residual
        if corrected is not None:
            corrected.append(sparse.csr_matrix(corrected_counts_tile(
                dense,
                library,
                intercept[indices],
                slope[indices],
                theta[indices],
                min_variance,
                scale_factor=cfg.corrected_scale_factor,
            )))
    timing["residual_statistics"] = time.perf_counter() - started

    started = time.perf_counter()
    hvg = np.argsort(-residual_variance, kind="stable")[:min(cfg.features, len(genes))]
    if full is None:
        residuals = np.empty((len(hvg), matrix.shape[1]), dtype=np.float64)
        for offset in range(0, len(hvg), cfg.tile_genes):
            indices = hvg[offset:offset + cfg.tile_genes]
            residuals[offset:offset + len(indices)] = cpu_tile(
                matrix[indices].toarray(),
                library,
                intercept[indices],
                slope[indices],
                theta[indices],
                min_variance,
                clip,
                dtype=np.float64,
            )[0]
        residual_rows = hvg
    else:
        residuals = full
        residual_rows = np.arange(len(genes))
    timing["selection_output"] = time.perf_counter() - started

    result = SCTResult(
        genes=genes,
        cells=cells,
        intercept=intercept,
        theta=theta,
        slope=slope,
        residual_mean=residual_mean,
        residual_variance=residual_variance,
        hvg_indices=hvg,
        residual_gene_indices=residual_rows,
        residuals=residuals,
        fit_gene_indices=eligible,
        fit_cell_indices=selected_cells,
        fit_intercept=fit_intercept,
        fit_theta=fit_theta,
        min_variance=min_variance,
        bandwidth=bandwidth,
        timing=timing,
        diagnostics={"fit_outliers": outliers},
    )
    if corrected is not None:
        result.corrected_umi = sparse.vstack(corrected, format="csr")
    return result
