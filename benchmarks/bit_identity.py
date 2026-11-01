from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from sctransform_cpu import SCTConfig, sctransform
from sctransform_cpu.r_compat_stats import row_statistics_r45


def _lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines()


def _float_bits(value) -> np.ndarray:
    array = np.ascontiguousarray(value, dtype="<f8")
    return array.view("<u8").reshape(-1)


def _compare_float(name: str, actual, oracle: Path) -> dict[str, object]:
    actual_array = np.ascontiguousarray(actual, dtype="<f8").reshape(-1)
    expected = np.fromfile(oracle / f"{name}.f64le", dtype="<f8")
    if actual_array.shape != expected.shape:
        return {
            "passed": False,
            "actual_elements": len(actual_array),
            "expected_elements": len(expected),
            "shape_mismatch": True,
        }
    actual_bits = _float_bits(actual_array)
    expected_bits = _float_bits(expected)
    mismatch = np.flatnonzero(actual_bits != expected_bits)
    finite = np.isfinite(actual_array) & np.isfinite(expected)
    delta = np.abs(actual_array[finite] - expected[finite])
    examples = []
    for index in mismatch[:5]:
        examples.append({
            "flat_index": int(index),
            "python": float(actual_array[index]),
            "r": float(expected[index]),
            "python_bits": f"0x{int(actual_bits[index]):016x}",
            "r_bits": f"0x{int(expected_bits[index]):016x}",
        })
    return {
        "elements": len(actual_array),
        "bit_mismatches": int(len(mismatch)),
        "bit_identical_fraction": float(1.0 - len(mismatch) / len(actual_array)),
        "max_abs_finite": float(delta.max()) if len(delta) else None,
        "first_mismatches": examples,
        "passed": len(mismatch) == 0,
    }


def _compare_i32(name: str, actual, oracle: Path) -> dict[str, object]:
    actual_array = np.ascontiguousarray(actual, dtype="<i4").reshape(-1)
    expected = np.fromfile(oracle / f"{name}.i32le", dtype="<i4")
    mismatch = np.flatnonzero(actual_array != expected) if actual_array.shape == expected.shape else np.arange(max(len(actual_array), len(expected)))
    return {
        "actual_elements": len(actual_array),
        "expected_elements": len(expected),
        "mismatches": int(len(mismatch)),
        "first_mismatch_indices": [int(value) for value in mismatch[:5]],
        "passed": actual_array.shape == expected.shape and len(mismatch) == 0,
    }


def compare_oracle(data: Path, oracle: Path) -> dict[str, object]:
    counts = sparse.load_npz(data / "counts.npz")
    genes = _lines(data / "gene_names.txt")
    cells = _lines(data / "cell_names.txt")
    result = sctransform(
        counts,
        gene_names=genes,
        cell_names=cells,
        config=SCTConfig(
            min_cells=5,
            fit_cells=None,
            fit_genes=None,
            features=len(genes),
            output="full",
            return_corrected_umi=True,
        ),
    )
    keep = np.asarray(counts.getnnz(axis=1)).ravel() >= 5
    retained = counts[keep].tocsr()
    means, variances, log_gmean = row_statistics_r45(retained)
    library = np.asarray(counts.sum(axis=0)).ravel()

    name_checks = {
        "genes": list(result.genes) == _lines(oracle / "genes.txt"),
        "cells": list(result.cells) == _lines(oracle / "cells.txt"),
        "fit_genes": list(result.genes[result.fit_gene_indices]) == _lines(oracle / "fit_genes.txt"),
        "fit_cells": list(result.cells[result.fit_cell_indices]) == _lines(oracle / "fit_cells.txt"),
        "hvg_order": list(result.genes[result.hvg_indices]) == _lines(oracle / "hvg_order.txt"),
    }
    floats = {
        "library": _compare_float("library", library, oracle),
        "amean": _compare_float("amean", means, oracle),
        "variance": _compare_float("variance", variances, oracle),
        "log_gmean": _compare_float("log_gmean", log_gmean, oracle),
        "fit_intercept": _compare_float("fit_intercept", result.fit_intercept, oracle),
        "fit_theta": _compare_float("fit_theta", result.fit_theta, oracle),
        "fit_slope": _compare_float("fit_slope", np.full(len(result.fit_theta), np.log(10.0)), oracle),
        "intercept": _compare_float("intercept", result.intercept, oracle),
        "theta": _compare_float("theta", result.theta, oracle),
        "slope": _compare_float("slope", np.full(len(result.theta), np.log(10.0)), oracle),
        "residual_mean": _compare_float("residual_mean", result.residual_mean, oracle),
        "residual_variance": _compare_float("residual_variance", result.residual_variance, oracle),
        "residual": _compare_float("residual", result.residuals, oracle),
    }
    integers = {
        "fit_outlier": _compare_i32("fit_outlier", result.diagnostics["fit_outliers"], oracle),
        "corrected_umi": _compare_i32("corrected_umi", result.corrected_umi.toarray(), oracle),
    }
    passed = (all(name_checks.values())
              and all(item["passed"] for item in floats.values())
              and all(item["passed"] for item in integers.values()))
    return {
        "contract": "exact IEEE-754 binary64 / int32 / ordered-title identity",
        "metadata": json.loads((oracle / "metadata.json").read_text(encoding="utf-8")),
        "name_checks": name_checks,
        "float_checks": floats,
        "integer_checks": integers,
        "passed": passed,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument("oracle", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--enforce", action="store_true")
    args = parser.parse_args()
    report = compare_oracle(args.data, args.oracle)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({
        "passed": report["passed"],
        "float_bit_mismatches": {
            name: value["bit_mismatches"]
            for name, value in report["float_checks"].items()
            if "bit_mismatches" in value
        },
        "integer_mismatches": {
            name: value["mismatches"]
            for name, value in report["integer_checks"].items()
        },
        "name_checks": report["name_checks"],
    }, indent=2))
    if args.enforce and not report["passed"]:
        raise SystemExit("bit-identity gate failed; inspect the JSON report")


if __name__ == "__main__":
    main()
