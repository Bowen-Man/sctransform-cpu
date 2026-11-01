import numpy as np
import pytest
from scipy import sparse

from sctransform_cpu import SCTConfig, SelectionManifest, sctransform
from sctransform_cpu.r_compat_stats import regularize_v2_r45
from sctransform_cpu.residuals import cpu_tile


def _counts(seed=19, genes=120, cells=80):
    rng = np.random.default_rng(seed)
    library = rng.lognormal(4.7, 0.3, size=cells)
    base = rng.lognormal(-4.2, 0.9, size=genes)
    mean = base[:, None] * library[None, :]
    return rng.negative_binomial(2.0, 2.0 / (2.0 + mean)).astype(np.int32)


def test_vst_contract_adapter_and_manifest_reuse():
    counts = _counts()
    genes = np.asarray([f"gene-{i}" for i in range(len(counts))])
    cells = np.asarray([f"cell-{i}" for i in range(counts.shape[1])])
    result = sctransform(
        sparse.csr_matrix(counts),
        gene_names=genes,
        cell_names=cells,
        config=SCTConfig(min_cells=2, output="full", features=40),
    )
    contract = result.to_vst_dict()
    assert list(contract) == [
        "y", "model_str", "model_pars", "model_pars_outliers",
        "model_pars_fit", "model_str_nonreg", "model_pars_nonreg",
        "arguments", "genes_log_gmean_step1", "cells_step1", "cell_attr",
        "gene_attr", "times",
    ]
    assert contract["y"]["values"].shape == (len(result.genes), len(cells))
    assert contract["model_pars_fit"]["values"].shape == (len(result.genes), 3)
    assert contract["gene_attr"]["residual_variance"].shape == (len(result.genes),)

    manifest = result.selection_manifest(42)
    repeated = sctransform(
        sparse.csr_matrix(counts),
        gene_names=genes,
        cell_names=cells,
        config=SCTConfig(min_cells=2, output="full", features=40),
        selection_manifest=manifest,
    )
    assert np.array_equal(repeated.fit_cell_indices, result.fit_cell_indices)
    assert np.array_equal(repeated.fit_gene_indices, result.fit_gene_indices)


def test_compact_result_refuses_full_r_contract():
    result = sctransform(
        sparse.csr_matrix(_counts()),
        config=SCTConfig(min_cells=2, output="hvg", features=20),
    )
    with pytest.raises(ValueError, match="output='full'"):
        result.to_vst_dict()


@pytest.mark.parametrize("bad", [
    np.asarray([[1.5, 0.0], [0.0, 2.0]]),
    np.asarray([[-1.0, 0.0], [0.0, 2.0]]),
    np.asarray([[np.inf, 0.0], [0.0, 2.0]]),
])
def test_invalid_raw_counts_are_rejected(bad):
    with pytest.raises(ValueError, match="raw UMI"):
        sctransform(bad, config=SCTConfig(min_cells=1))


def test_zero_library_and_all_filtered_inputs_are_rejected():
    with pytest.raises(ValueError, match="zero-library"):
        sctransform(np.asarray([[1, 0], [2, 0]]), config=SCTConfig(min_cells=1))
    with pytest.raises(ValueError, match="too few retained"):
        sctransform(np.eye(4), config=SCTConfig(min_cells=2))


def test_poisson_boundaries_and_residual_clipping():
    fit_count = 30
    fit_x = np.linspace(-2.0, 1.0, fit_count)
    log_gmean = np.concatenate((fit_x, [-3.0, -4.0]))
    means = np.concatenate((np.linspace(0.2, 5.0, fit_count), [1.0, 5e-4]))
    variances = np.concatenate((means[:fit_count] + 0.5, [1.0, 2e-3]))
    intercept = -2.0 + 0.7 * fit_x
    theta = np.linspace(1.5, 4.0, fit_count)
    fitted_intercept, fitted_theta, fitted_slope, _, _ = regularize_v2_r45(
        log_gmean,
        np.arange(fit_count),
        intercept,
        theta,
        means,
        variances,
        mean_filtered_library=100.0,
    )
    assert np.isinf(fitted_theta[-2:]).all()
    assert fitted_intercept[-2] == np.log(means[-2]) - np.log(100.0)
    assert fitted_intercept[-1] == np.log(means[-1]) - np.log(100.0)
    assert np.all(fitted_slope[-2:] == np.log(10.0))

    counts = np.asarray([[0.0, 0.0, 0.0, 1e6]])
    residual, _, _ = cpu_tile(
        counts,
        np.asarray([10.0, 10.0, 10.0, 10.0]),
        np.asarray([-2.0]),
        np.asarray([np.log(10.0)]),
        np.asarray([2.0]),
        min_variance=0.01,
        clip=2.0,
        dtype=np.float64,
    )
    assert np.max(residual) == 2.0
    assert np.min(residual) >= -2.0
