from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import spearmanr

from sctransform_cpu import SCTConfig, sctransform


TOLERANCE = {
    "intercept": (1e-7, 1e-7),
    "theta": (1e-8, 1e-5),
    "residual": (1e-6, 1e-5),
    "residual_variance": (1e-7, 1e-5),
    "slope": (1e-13, 0.0),
}


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def values(rows, column):
    return np.asarray([float(row[column]) for row in rows])


def compare(actual, expected, key):
    actual, expected = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    atol, rtol = TOLERANCE[key]
    finite = np.isfinite(actual) & np.isfinite(expected)
    nonfinite_mismatch = int(np.count_nonzero(np.isfinite(actual) != np.isfinite(expected)))
    delta = np.abs(actual[finite] - expected[finite])
    allowed = atol + rtol * np.abs(expected[finite])
    return {
        "atol": atol,
        "rtol": rtol,
        "finite_pairs": int(finite.sum()),
        "nonfinite_mismatch": nonfinite_mismatch,
        "tolerance_failures": int(np.count_nonzero(delta > allowed)) + nonfinite_mismatch,
        "mae": float(delta.mean()) if len(delta) else None,
        "rmse": float(np.sqrt(np.mean(delta * delta))) if len(delta) else None,
        "max_abs": float(delta.max()) if len(delta) else None,
        "max_scaled_error": float(np.max(delta / allowed)) if len(delta) else None,
        "passed": nonfinite_mismatch == 0 and bool(np.all(delta <= allowed)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument("oracle", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    counts = sparse.load_npz(args.data / "counts.npz")
    genes = (args.data / "gene_names.txt").read_text(encoding="utf-8").splitlines()
    cells = (args.data / "cell_names.txt").read_text(encoding="utf-8").splitlines()
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
    fit = read_csv(args.oracle / "fit.csv")
    all_genes = read_csv(args.oracle / "genes.csv")
    residual_rows = read_csv(args.oracle / "residuals.csv")
    corrected_rows = read_csv(args.oracle / "corrected_umi.csv")
    if [row["gene"] for row in all_genes] != list(result.genes):
        raise ValueError("retained gene order differs")
    if [row["gene"] for row in fit] != list(result.genes[result.fit_gene_indices]):
        raise ValueError("fit gene order differs")
    r_residual = np.asarray([[float(row[cell]) for cell in cells] for row in residual_rows])
    r_corrected = np.asarray([[float(row[cell]) for cell in cells] for row in corrected_rows])
    r_variance = values(all_genes, "residual_variance")
    r_outlier = np.asarray([row["outlier"] == "TRUE" for row in fit])
    checks = {
        "step1_intercept": compare(result.fit_intercept, values(fit, "(Intercept)"), "intercept"),
        "step1_theta": compare(result.fit_theta, values(fit, "theta"), "theta"),
        "step1_slope": compare(np.full(len(fit), np.log(10.0)), values(fit, "log_umi"), "slope"),
        "regularized_intercept": compare(result.intercept, values(all_genes, "(Intercept)"), "intercept"),
        "regularized_theta": compare(result.theta, values(all_genes, "theta"), "theta"),
        "regularized_slope": compare(result.slope, values(all_genes, "log_umi"), "slope"),
        "residual": compare(result.residuals, r_residual, "residual"),
        "residual_variance": compare(result.residual_variance, r_variance, "residual_variance"),
    }
    corrected_equal = bool(np.array_equal(result.corrected_umi.toarray(), r_corrected))
    r_order = np.argsort(-r_variance, kind="stable")
    hvg = {}
    for size in (100, 200, 300):
        left = set(result.genes[result.hvg_indices[:size]])
        right = set(result.genes[r_order[:size]])
        hvg[str(size)] = {
            "intersection": len(left & right),
            "jaccard": len(left & right) / len(left | right),
            "passed": left == right,
        }
    outlier_mismatch = int(np.count_nonzero(result.diagnostics["fit_outliers"] != r_outlier))
    passed = (all(item["passed"] for item in checks.values())
              and corrected_equal and outlier_mismatch == 0
              and all(item["passed"] for item in hvg.values()))
    payload = {
        "scope": "PBMC3K-derived 800 gene x 240 cell direct unbatched v2 path",
        "reference": "R 4.5.2, sctransform 0.4.3, glmGamPoi 1.22.0, matrixStats 1.5.0",
        "retained_genes": len(result.genes),
        "fitted_genes": len(result.fit_gene_indices),
        "fitted_cells": len(result.fit_cell_indices),
        "checks": checks,
        "outlier_mismatch": outlier_mismatch,
        "corrected_umi_exact": corrected_equal,
        "residual_variance_spearman": float(spearmanr(result.residual_variance, r_variance).statistic),
        "hvg": hvg,
        "python_stage_seconds_single_quality_run": result.timing,
        "passed": passed,
        "not_proven": [
            "other inputs and numerical boundaries",
            "implicit cross-language sampling",
            "batch/covariates and complete Seurat workflow",
            "other operating systems and package versions",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    if not passed:
        raise SystemExit("R pointwise quality gate failed")


if __name__ == "__main__":
    main()
