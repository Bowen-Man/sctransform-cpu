from __future__ import annotations

import os
from pathlib import Path
import subprocess

import numpy as np
import pytest
from scipy import sparse

from sctransform_cpu import SCTConfig, sctransform


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.skipif(
    not os.environ.get("SCT_TEST_RSCRIPT"),
    reason="set SCT_TEST_RSCRIPT to run the live R design gate",
)
@pytest.mark.parametrize(
    "mode", ["covariate", "batch", "covariate_corrected", "batch_nonreg"])
def test_live_r_design_bit_identity(tmp_path, mode):
    oracle = tmp_path / mode
    subprocess.run([
        os.environ["SCT_TEST_RSCRIPT"],
        str(ROOT / "benchmarks" / "design_reference.R"),
        str(ROOT / "benchmark_data"), str(oracle), mode,
    ], check=True)

    counts = sparse.load_npz(ROOT / "benchmark_data" / "counts.npz")
    genes = np.loadtxt(ROOT / "benchmark_data" / "genes.tsv", dtype=str, skiprows=1)
    cells = np.loadtxt(ROOT / "benchmark_data" / "cells.tsv", dtype=str, skiprows=1)
    attributes = {
        "percent": (np.arange(counts.shape[1]) % 17) / 17,
        "batch": np.resize(np.asarray(["A", "B"]), counts.shape[1]),
    }
    is_batch = mode in {"batch", "batch_nonreg"}
    result = sctransform(
        counts, gene_names=genes, cell_names=cells, cell_attributes=attributes,
        config=SCTConfig(
            min_cells=5, fit_cells=None, fit_genes=None, output="full",
            latent_variables=("log_umi", "percent") if not is_batch
            else ("log_umi",),
            latent_variables_nonreg=("percent",) if mode == "batch_nonreg"
            else (),
            batch_variable="batch" if is_batch else None,
            return_corrected_umi=mode in {
                "covariate_corrected", "batch_nonreg"},
            corrected_scale_factor=1234.0
            if mode == "covariate_corrected" else None,
        ),
    )

    fit_genes = (oracle / "fit_genes.txt").read_text().splitlines()
    out_genes = (oracle / "genes.txt").read_text().splitlines()
    columns = (oracle / "columns.txt").read_text().splitlines()
    assert result.genes.tolist() == out_genes
    assert result.genes[result.fit_gene_indices].tolist() == fit_genes
    assert ["theta", *result.coefficient_names] == columns

    expected = {
        "model_pars": np.fromfile(oracle / "model_pars.f64le", "<f8").reshape(
            len(fit_genes), len(columns)),
        "model_pars_fit": np.fromfile(
            oracle / "model_pars_fit.f64le", "<f8").reshape(len(out_genes), len(columns)),
        "residual": np.fromfile(oracle / "residual.f64le", "<f8").reshape(
            len(out_genes), len(cells)),
        "residual_mean": np.fromfile(oracle / "residual_mean.f64le", "<f8"),
        "residual_variance": np.fromfile(oracle / "residual_variance.f64le", "<f8"),
        "outliers": np.fromfile(oracle / "outliers.f64le", "<f8"),
    }
    actual = {
        "model_pars": np.column_stack((result.fit_theta, result.fit_coefficients)),
        "model_pars_fit": np.column_stack((result.theta, result.coefficients)),
        "residual": result.residuals,
        "residual_mean": result.residual_mean,
        "residual_variance": result.residual_variance,
        "outliers": result.diagnostics["fit_outliers"].astype(float),
    }
    for name in expected:
        assert np.array_equal(
            actual[name].view(np.uint64), expected[name].view(np.uint64)), name

    nonreg_columns = oracle / "nonreg_columns.txt"
    if nonreg_columns.exists():
        names = nonreg_columns.read_text().splitlines()
        values = np.fromfile(
            oracle / "model_pars_nonreg.f64le", "<f8").reshape(
                len(out_genes), len(names))
        assert result.nonreg_coefficient_names == tuple(names)
        assert np.array_equal(
            result.nonreg_coefficients.view(np.uint64), values.view(np.uint64))
    else:
        assert result.nonreg_coefficients is None

    corrected = oracle / "corrected_umi.i32le"
    if corrected.exists():
        values = np.fromfile(corrected, "<i4").reshape(len(out_genes), len(cells))
        assert np.array_equal(result.corrected_umi.toarray(), values)
    else:
        assert result.corrected_umi is None
