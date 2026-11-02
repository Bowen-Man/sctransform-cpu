"""Minimal CPU-only implementation of the validated SCT v2 paths."""
from __future__ import annotations

from dataclasses import dataclass, field
import time

import numpy as np
from scipy import sparse

from .corrected_counts import corrected_counts_design_tile, corrected_counts_tile
from .design import build_design, correction_design
from .fit import fit_batch, fit_general_batch
from .r_compat_stats import (
    regularize_design_v2_r45,
    regularize_v2_r45,
    row_statistics_r45,
)
from .residuals import cpu_design_tile, cpu_tile
from .port_native import general_linear_native, nonreg_fit_native
from .sampling import (
    RMersenneTwister,
    SelectionManifest,
    r_nrd_density_weights,
)


@dataclass(frozen=True)
class SCTConfig:
    min_cells: int = 5
    fit_cells: int | None = 2000
    fit_genes: int | None = 2000
    random_seed: int | None = None
    features: int = 3000
    output: str = "hvg"
    tile_genes: int = 256
    bandwidth_adjust: float = 3.0
    return_corrected_umi: bool = False
    corrected_scale_factor: float | None = None
    latent_variables: tuple[str, ...] = ("log_umi",)
    latent_variables_nonreg: tuple[str, ...] = ()
    batch_variable: str | None = None


@dataclass
class SCTResult:
    genes: np.ndarray
    cells: np.ndarray
    intercept: np.ndarray | None
    theta: np.ndarray
    slope: np.ndarray | None
    residual_mean: np.ndarray
    residual_variance: np.ndarray
    hvg_indices: np.ndarray
    residual_gene_indices: np.ndarray
    residuals: np.ndarray
    fit_gene_indices: np.ndarray
    fit_cell_indices: np.ndarray
    fit_intercept: np.ndarray | None
    fit_theta: np.ndarray
    library_size: np.ndarray
    log_gmean: np.ndarray
    amean: np.ndarray
    variance: np.ndarray
    detection_rate: np.ndarray
    min_variance: float
    bandwidth: float
    timing: dict[str, float] = field(default_factory=dict)
    diagnostics: dict[str, object] = field(default_factory=dict)
    corrected_umi: sparse.csr_matrix | None = None
    coefficient_names: tuple[str, ...] = ("(Intercept)", "log_umi")
    coefficients: np.ndarray | None = None
    fit_coefficients: np.ndarray | None = None
    model_str: str = "y ~ log_umi"
    nonreg_coefficient_names: tuple[str, ...] = ()
    nonreg_coefficients: np.ndarray | None = None
    model_str_nonreg: str = ""
    semantics: str = (
        "platform-pinned Windows x64 or Linux x86_64 R 4.5.2 / "
        "sctransform 0.4.3 / "
        "glmGamPoi 1.22.0 / matrixStats 1.5.0; regularized numeric and "
        "single categorical-batch v2 designs, finite-theta nonreg and "
        "corrected UMI"
    )

    def selection_manifest(self, seed: int) -> SelectionManifest:
        """Return ordered fit titles that can be reused across implementations."""
        return SelectionManifest(
            int(seed),
            tuple(self.cells[self.fit_cell_indices]),
            tuple(self.genes[self.fit_gene_indices]),
        )

    @staticmethod
    def _named_matrix(rows, columns, values) -> dict[str, np.ndarray]:
        return {
            "row_names": np.asarray(rows, dtype=str),
            "column_names": np.asarray(columns, dtype=str),
            "values": np.asarray(values),
        }

    def to_vst_dict(self) -> dict[str, object]:
        """Expose the targeted R ``vst`` return contract with Python containers.

        R matrices are represented by ``row_names``, ``column_names`` and
        ``values``.  A full R-like ``y`` requires ``output='full'``; compact
        HVG output remains available through the dataclass itself.
        """
        if len(self.residual_gene_indices) != len(self.genes):
            raise ValueError("to_vst_dict requires SCTConfig(output='full')")
        fitted_coefficients = self.fit_coefficients
        regularized_coefficients = self.coefficients
        if fitted_coefficients is None or regularized_coefficients is None:
            fitted_coefficients = np.column_stack(
                (self.fit_intercept, np.full(len(self.fit_theta), np.log(10.0))))
            regularized_coefficients = np.column_stack((self.intercept, self.slope))
        fitted = np.column_stack((self.fit_theta, fitted_coefficients))
        regularized = np.column_stack((self.theta, regularized_coefficients))
        columns = ("theta", *self.coefficient_names)
        result: dict[str, object] = {
            "y": self._named_matrix(self.genes, self.cells, self.residuals),
            "model_str": self.model_str,
            "model_pars": self._named_matrix(
                self.genes[self.fit_gene_indices],
                columns,
                fitted,
            ),
            "model_pars_outliers": np.asarray(
                self.diagnostics["fit_outliers"], dtype=bool),
            "model_pars_fit": self._named_matrix(
                self.genes,
                columns,
                regularized,
            ),
            "model_str_nonreg": self.model_str_nonreg,
            "model_pars_nonreg": (
                self._named_matrix(
                    self.genes, self.nonreg_coefficient_names,
                    self.nonreg_coefficients)
                if self.nonreg_coefficients is not None else None),
            "arguments": dict(self.diagnostics["arguments"]),
            "genes_log_gmean_step1": self.log_gmean[self.fit_gene_indices].copy(),
            "cells_step1": self.cells[self.fit_cell_indices].copy(),
            "cell_attr": None,
            "gene_attr": {
                "row_names": self.genes.copy(),
                "detection_rate": self.detection_rate.copy(),
                "gmean": np.power(10.0, self.log_gmean),
                "amean": self.amean.copy(),
                "variance": self.variance.copy(),
                "residual_mean": self.residual_mean.copy(),
                "residual_variance": self.residual_variance.copy(),
            },
            "times": dict(self.timing),
        }
        if self.corrected_umi is not None:
            result["umi_corrected"] = self.corrected_umi
        return result


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
                fit_cell_indices=None, fit_gene_names=None,
                selection_manifest: SelectionManifest | None = None,
                cell_attributes=None) -> SCTResult:
    """Transform raw UMI counts in gene x cell orientation.

    Exact comparison can either reuse an ordered ``SelectionManifest`` or use
    ``random_seed`` to reproduce R 4.5's default sampling stream.
    """
    cfg = config or SCTConfig()
    if cfg.min_cells < 1 or cfg.features < 1 or cfg.tile_genes < 1:
        raise ValueError("min_cells, features and tile_genes must be positive")
    if cfg.bandwidth_adjust <= 0 or cfg.output not in {"hvg", "full"}:
        raise ValueError("invalid bandwidth_adjust or output")
    if any(value is not None and value < 1 for value in (cfg.fit_cells, cfg.fit_genes)):
        raise ValueError("fit_cells and fit_genes must be positive or None")
    general_design = (cfg.latent_variables != ("log_umi",)
                      or cfg.batch_variable is not None)
    has_nonreg = bool(cfg.latent_variables_nonreg)

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
    input_matrix = matrix
    library = np.asarray(matrix.sum(axis=0)).ravel()
    keep = np.asarray(matrix.getnnz(axis=1)).ravel() >= cfg.min_cells
    matrix, genes = matrix[keep].tocsr(), genes[keep]
    if matrix.shape[0] < 2 or np.any(library <= 0):
        raise ValueError("too few retained genes or a zero-library cell")
    retained_mean_library = float(np.asarray(matrix.sum(axis=0)).mean())
    means, variances, log_gmean = row_statistics_r45(matrix)
    detection_rate = np.asarray(matrix.getnnz(axis=1)).ravel() / matrix.shape[1]
    timing = {"input_statistics": time.perf_counter() - started}

    started = time.perf_counter()
    if selection_manifest is not None:
        if fit_cell_indices is not None or fit_gene_names is not None:
            raise ValueError("selection_manifest cannot be combined with explicit selections")
        cell_lookup = {name: index for index, name in enumerate(cells)}
        if any(name not in cell_lookup for name in selection_manifest.cells):
            raise ValueError("selection manifest contains an unknown cell")
        fit_cell_indices = np.asarray(
            [cell_lookup[name] for name in selection_manifest.cells], dtype=int)
        fit_gene_names = np.asarray(selection_manifest.genes, dtype=str)

    sampler = None
    if fit_cell_indices is None:
        requested_cells = matrix.shape[1] if cfg.fit_cells is None else cfg.fit_cells
        if requested_cells < matrix.shape[1]:
            if cfg.random_seed is None:
                raise ValueError("implicit sampling requires SCTConfig(random_seed=...)")
            sampler = RMersenneTwister(cfg.random_seed)
            selected_cells = sampler.sample_uniform(matrix.shape[1], requested_cells)
        else:
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
        if cfg.random_seed is None:
            raise ValueError("implicit sampling requires SCTConfig(random_seed=...)")
        if sampler is None:
            sampler = RMersenneTwister(cfg.random_seed)
        selected = sampler.sample_weighted(
            r_nrd_density_weights(log_gmean[eligible]), cfg.fit_genes)
        eligible = eligible[selected]
    if len(eligible) < 2:
        raise ValueError("fewer than two overdispersed genes are available for fitting")
    timing["selection"] = time.perf_counter() - started

    started = time.perf_counter()
    design = build_design(
        input_matrix, library, cfg.latent_variables, cfg.batch_variable,
        cfg.latent_variables_nonreg, cell_attributes, selected_cells)
    fit_theta = np.empty(len(eligible))
    fit_coefficients = np.empty((len(eligible), len(design.coefficient_names)))
    if not general_design:
        fit_intercept = fit_coefficients[:, 0]
        for offset in range(0, len(eligible), 500):
            indices = eligible[offset:offset + 500]
            fit_intercept[offset:offset + len(indices)], fit_theta[offset:offset + len(indices)] = fit_batch(
                sampled[indices].toarray(), library[selected_cells],
                means[indices], variances[indices])
        fit_coefficients[:, 1] = np.log(10.0)
    else:
        for offset in range(0, len(eligible), 500):
            indices = eligible[offset:offset + 500]
            raw_coefficients, values = fit_general_batch(
                sampled[indices].toarray(), design.fit, design.offset,
                means[indices], variances[indices])
            fit_theta[offset:offset + len(indices)] = values
            if cfg.batch_variable is None:
                fit_coefficients[offset:offset + len(indices), 0] = raw_coefficients[:, 0]
                fit_coefficients[offset:offset + len(indices), 1] = np.log(10.0)
                fit_coefficients[offset:offset + len(indices), 2:] = raw_coefficients[:, 1:]
            else:
                fit_coefficients[offset:offset + len(indices)] = raw_coefficients
        fit_theta[np.isfinite(fit_theta) & (fit_theta < 1e-7)] = 1e-7
        fit_intercept = (fit_coefficients[:, design.coefficient_names.index("(Intercept)")]
                         if "(Intercept)" in design.coefficient_names else None)
    timing["fit"] = time.perf_counter() - started

    started = time.perf_counter()
    if not general_design:
        intercept, theta, slope, bandwidth, outliers = regularize_v2_r45(
            log_gmean, eligible, fit_intercept, fit_theta, means, variances,
            retained_mean_library, cfg.bandwidth_adjust)
        coefficients = np.column_stack((intercept, slope))
    else:
        coefficients, theta, bandwidth, outliers = regularize_design_v2_r45(
            log_gmean, eligible, fit_coefficients, fit_theta,
            design.coefficient_names, means, variances, retained_mean_library,
            cfg.bandwidth_adjust, counts=matrix,
            batch=design.attributes.get(cfg.batch_variable),
            fit_cells=selected_cells, batch_variable=cfg.batch_variable)
        intercept = (coefficients[:, design.coefficient_names.index("(Intercept)")]
                     if "(Intercept)" in design.coefficient_names else None)
        slope = (coefficients[:, design.coefficient_names.index("log_umi")]
                 if "log_umi" in design.coefficient_names else None)
    timing["regularization"] = time.perf_counter() - started

    started = time.perf_counter()
    nonreg_coefficients = None
    final_design = design.full
    final_coefficients = coefficients
    if design.nonreg_full is not None:
        nonreg_coefficients = np.empty(
            (len(genes), len(design.nonreg_coefficient_names)),
            dtype=np.float64)
        for offset in range(0, len(genes), 500):
            indices = np.arange(offset, min(offset + 500, len(genes)))
            base_linear = general_linear_native(
                coefficients[indices], design.full)
            fitted_nonreg = nonreg_fit_native(
                matrix[indices].toarray(), design.nonreg_full,
                base_linear, theta[indices])
            nonreg_coefficients[indices] = fitted_nonreg
        final_design = np.column_stack((design.full, design.nonreg_full))
        final_coefficients = np.column_stack(
            (coefficients, nonreg_coefficients))
    timing["nonregularized_fit"] = time.perf_counter() - started

    started = time.perf_counter()
    clip = float(np.sqrt(matrix.shape[1]))
    min_variance = float((np.median(matrix.data) / 5.0) ** 2)
    residual_mean = np.empty(len(genes))
    residual_variance = np.empty(len(genes))
    full = np.empty(matrix.shape, dtype=np.float64) if cfg.output == "full" else None
    corrected = [] if cfg.return_corrected_umi else None
    target_design = (correction_design(
        design, cfg.latent_variables, cfg.batch_variable,
        cfg.corrected_scale_factor) if corrected is not None else None)
    for offset in range(0, len(genes), cfg.tile_genes):
        indices = np.arange(offset, min(offset + cfg.tile_genes, len(genes)))
        dense = matrix[indices].toarray()
        if general_design or has_nonreg:
            residual, mean, variance = cpu_design_tile(
                dense, final_design, final_coefficients[indices], theta[indices],
                min_variance, clip)
        else:
            residual, mean, variance = cpu_tile(
                dense, library, intercept[indices], slope[indices], theta[indices],
                min_variance, clip, dtype=np.float64)
        residual_mean[indices], residual_variance[indices] = mean, variance
        if full is not None:
            full[indices] = residual
        if corrected is not None:
            if general_design or has_nonreg:
                corrected_values = corrected_counts_design_tile(
                    dense, final_design, final_coefficients[indices],
                    target_design, coefficients[indices], theta[indices],
                    min_variance)
            else:
                corrected_values = corrected_counts_tile(
                    dense, library, intercept[indices], slope[indices],
                    theta[indices], min_variance,
                    scale_factor=cfg.corrected_scale_factor)
            corrected.append(sparse.csr_matrix(corrected_values))
    timing["residual_statistics"] = time.perf_counter() - started

    started = time.perf_counter()
    hvg = np.argsort(-residual_variance, kind="stable")[:min(cfg.features, len(genes))]
    if full is None:
        residuals = np.empty((len(hvg), matrix.shape[1]), dtype=np.float64)
        for offset in range(0, len(hvg), cfg.tile_genes):
            indices = hvg[offset:offset + cfg.tile_genes]
            dense = matrix[indices].toarray()
            if general_design or has_nonreg:
                residuals[offset:offset + len(indices)] = cpu_design_tile(
                    dense, final_design, final_coefficients[indices], theta[indices],
                    min_variance, clip)[0]
            else:
                residuals[offset:offset + len(indices)] = cpu_tile(
                    dense, library, intercept[indices], slope[indices],
                    theta[indices], min_variance, clip, dtype=np.float64)[0]
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
        library_size=library,
        log_gmean=log_gmean,
        amean=means,
        variance=variances,
        detection_rate=detection_rate,
        min_variance=min_variance,
        bandwidth=bandwidth,
        coefficient_names=design.coefficient_names,
        coefficients=coefficients,
        fit_coefficients=fit_coefficients,
        model_str=design.model_str,
        nonreg_coefficient_names=design.nonreg_coefficient_names,
        nonreg_coefficients=nonreg_coefficients,
        model_str_nonreg=design.model_str_nonreg,
        timing=timing,
        diagnostics={
            "fit_outliers": outliers,
            "arguments": {
                "vst.flavor": "v2",
                "method": "glmGamPoi_offset",
                "n_cells": cfg.fit_cells,
                "n_genes": cfg.fit_genes,
                "min_cells": cfg.min_cells,
                "residual_type": "pearson",
                "min_variance": "umi_median",
                "set_min_var": min_variance,
                "bw_adjust": cfg.bandwidth_adjust,
                "exclude_poisson": True,
                "return_corrected_umi": cfg.return_corrected_umi,
                "latent_var": cfg.latent_variables,
                "latent_var_nonreg": cfg.latent_variables_nonreg,
                "batch_var": cfg.batch_variable,
            },
        },
    )
    if corrected is not None:
        result.corrected_umi = sparse.vstack(corrected, format="csr")
    return result
