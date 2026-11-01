from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy import sparse

from sctransform_cpu import SCTConfig, sctransform


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    counts = sparse.load_npz(args.data / "counts.npz")
    genes = (args.data / "gene_names.txt").read_text(encoding="utf-8").splitlines()
    cells = (args.data / "cell_names.txt").read_text(encoding="utf-8").splitlines()
    wall_start = time.perf_counter()
    cpu_start = time.process_time()
    result = sctransform(
        counts,
        gene_names=genes,
        cell_names=cells,
        config=SCTConfig(
            min_cells=5,
            fit_cells=None,
            fit_genes=None,
            features=300,
            output="full",
            return_corrected_umi=True,
        ),
    )
    payload = {
        "implementation": "sctransform-cpu",
        "compute_wall_seconds": time.perf_counter() - wall_start,
        "compute_cpu_seconds": time.process_time() - cpu_start,
        "stage_seconds": result.timing,
        "retained_genes": len(result.genes),
        "fitted_genes": len(result.fit_gene_indices),
        "fitted_cells": len(result.fit_cell_indices),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
