from __future__ import annotations

import os
from pathlib import Path
import subprocess

import numpy as np
import pytest
from scipy import sparse
from scipy.io import mmwrite

from benchmarks.bit_identity import compare_oracle


ROOT = Path(__file__).resolve().parents[1]


def _write_edge_fixture(path: Path) -> None:
    rng = np.random.default_rng(20261102)
    cells = 96
    ordinary = 176
    depth = rng.lognormal(0.0, 0.4, size=cells)
    abundance = rng.lognormal(1.5, 0.7, size=ordinary)
    mean = abundance[:, None] * depth[None, :]
    base = rng.negative_binomial(1.25, 1.25 / (1.25 + mean)).astype(np.int32)

    exact_boundary = np.zeros(cells, dtype=np.int32)
    exact_boundary[:77] = 5  # sample variance equals mean for n=96
    rare_kept = np.zeros(cells, dtype=np.int32)
    rare_kept[[1, 7, 19, 43, 89]] = 1
    rare_filtered = np.zeros(cells, dtype=np.int32)
    rare_filtered[[2, 11, 31, 71]] = 1
    extreme = np.ones(cells, dtype=np.int32)
    extreme[0] = 1_000_000
    counts = np.vstack((
        base,
        np.ones(cells, dtype=np.int32),
        np.full(cells, 2, dtype=np.int32),
        exact_boundary,
        rare_kept,
        rare_filtered,
        extreme,
    ))

    path.mkdir()
    sparse.save_npz(path / "counts.npz", sparse.csr_matrix(counts))
    mmwrite(path / "counts.mtx", sparse.coo_matrix(counts))
    genes = [f"gene-{index:04d}" for index in range(len(counts))]
    cells_names = [f"cell-{index:03d}" for index in range(cells)]
    (path / "genes.tsv").write_text(
        "gene\n" + "\n".join(genes) + "\n", encoding="utf-8")
    (path / "cells.tsv").write_text(
        "cell\n" + "\n".join(cells_names) + "\n", encoding="utf-8")
    (path / "gene_names.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
    (path / "cell_names.txt").write_text("\n".join(cells_names) + "\n", encoding="utf-8")


@pytest.mark.skipif(
    not os.environ.get("SCT_TEST_RSCRIPT"),
    reason="set SCT_TEST_RSCRIPT to run the live R adversarial gate",
)
def test_live_r_edge_fixture_bit_identity(tmp_path):
    data = tmp_path / "edge-data"
    oracle = tmp_path / "edge-oracle"
    _write_edge_fixture(data)
    subprocess.run([
        os.environ["SCT_TEST_RSCRIPT"],
        str(ROOT / "benchmarks" / "bit_reference.R"),
        str(data),
        str(oracle),
    ], check=True)
    report = compare_oracle(data, oracle)
    assert report["passed"], report
