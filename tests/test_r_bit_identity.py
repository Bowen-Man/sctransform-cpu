from __future__ import annotations

import os
from pathlib import Path
import subprocess

import json
import numpy as np
import pytest
from scipy import sparse
from scipy.io import mmwrite

from benchmarks.bit_identity import compare_oracle
from sctransform_cpu import SCTConfig, make_selection_manifest, sctransform


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.skipif(
    not os.environ.get("SCT_TEST_RSCRIPT"),
    reason="set SCT_TEST_RSCRIPT to run the live R bit-identity gate",
)
def test_live_r_bit_identity(tmp_path):
    rscript = os.environ["SCT_TEST_RSCRIPT"]
    oracle = tmp_path / "oracle"
    subprocess.run([
        rscript,
        str(ROOT / "benchmarks" / "bit_reference.R"),
        str(ROOT / "benchmark_data"),
        str(oracle),
    ], check=True)
    report = compare_oracle(ROOT / "benchmark_data", oracle)
    assert report["passed"], report


@pytest.mark.skipif(
    not os.environ.get("SCT_TEST_RSCRIPT"),
    reason="set SCT_TEST_RSCRIPT to run the live R sampling gate",
)
def test_live_r_implicit_selection_identity(tmp_path):
    rng = np.random.default_rng(20261102)
    cells = 60
    genes = 240
    depth = rng.lognormal(0.0, 0.25, size=cells)
    abundance = np.linspace(2.0, 18.0, genes)
    mean = abundance[:, None] * depth[None, :]
    counts = rng.negative_binomial(1.5, 1.5 / (1.5 + mean)).astype(np.int32)
    gene_names = np.asarray([f"gene-{index:04d}" for index in range(genes)])
    cell_names = np.asarray([f"cell-{index:03d}" for index in range(cells)])

    data = tmp_path / "data"
    data.mkdir()
    mmwrite(data / "counts.mtx", sparse.coo_matrix(counts))
    (data / "genes.tsv").write_text(
        "gene\n" + "\n".join(gene_names) + "\n", encoding="utf-8")
    (data / "cells.tsv").write_text(
        "cell\n" + "\n".join(cell_names) + "\n", encoding="utf-8")

    reference_path = tmp_path / "selection.json"
    subprocess.run([
        os.environ["SCT_TEST_RSCRIPT"],
        str(ROOT / "benchmarks" / "selection_reference.R"),
        str(data),
        str(reference_path),
        "42",
        "40",
        "40",
    ], check=True)
    reference = json.loads(reference_path.read_text(encoding="utf-8"))

    manifest = make_selection_manifest(
        sparse.csr_matrix(counts),
        gene_names=gene_names,
        cell_names=cell_names,
        seed=42,
        n_cells=40,
        n_genes=40,
        min_cells=5,
    )
    assert list(manifest.cells) == reference["cells"]
    assert list(manifest.genes) == reference["genes"]

    result = sctransform(
        sparse.csr_matrix(counts),
        gene_names=gene_names,
        cell_names=cell_names,
        config=SCTConfig(
            min_cells=5,
            fit_cells=40,
            fit_genes=40,
            random_seed=42,
            output="hvg",
            features=20,
        ),
    )
    assert result.cells[result.fit_cell_indices].tolist() == reference["cells"]
    assert result.genes[result.fit_gene_indices].tolist() == reference["genes"]
