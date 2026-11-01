"""Deterministic CPU profile for the fixed SCT regression workload.

The input is loaded and one warm-up is completed before profiling. This keeps
import, file I/O and one-time DLL loading out of the algorithm profile.
"""
from __future__ import annotations

import argparse
import cProfile
import json
from pathlib import Path
import pstats

from scipy import sparse

from sctransform_cpu import SCTConfig, sctransform


def run(counts, genes, cells):
    return sctransform(
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--limit", type=int, default=25)
    args = parser.parse_args()

    counts = sparse.load_npz(args.data / "counts.npz")
    genes = (args.data / "gene_names.txt").read_text(
        encoding="utf-8").splitlines()
    cells = (args.data / "cell_names.txt").read_text(
        encoding="utf-8").splitlines()

    run(counts, genes, cells)  # warm-up: imports, allocations and DLL load
    profiler = cProfile.Profile()
    profiler.enable()
    result = run(counts, genes, cells)
    profiler.disable()

    stats = pstats.Stats(profiler)
    project_rows = []
    for (filename, line, function), (primitive, calls, self_time,
                                      cumulative, _callers) in stats.stats.items():
        normalized = filename.replace("\\", "/")
        if "/src/sctransform_cpu/" not in normalized:
            continue
        project_rows.append({
            "file": Path(filename).name,
            "line": line,
            "function": function,
            "primitive_calls": primitive,
            "total_calls": calls,
            "self_seconds": self_time,
            "cumulative_seconds": cumulative,
        })
    project_rows.sort(key=lambda item: item["cumulative_seconds"], reverse=True)

    payload = {
        "method": "one unprofiled warm-up, then cProfile around sctransform only",
        "warning": (
            "cProfile perturbs runtime; cumulative rows are nested and must not be added"
        ),
        "profiled_total_calls": stats.total_calls,
        "profiled_primitive_calls": stats.prim_calls,
        "profiled_total_seconds": stats.total_tt,
        "stage_seconds_under_profile": result.timing,
        "top_project_functions": project_rows[:args.limit],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
