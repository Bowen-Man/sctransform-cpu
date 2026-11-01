from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import tempfile
import time

import psutil


ROOT = Path(__file__).resolve().parents[1]


def monitored(command, output: Path, env):
    started = time.perf_counter()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, env=env)
    peak = 0
    while process.poll() is None:
        try:
            root = psutil.Process(process.pid)
            members = [root, *root.children(recursive=True)]
            peak = max(peak, sum(member.memory_info().rss for member in members
                                 if member.is_running()))
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        time.sleep(0.01)
    stdout, stderr = process.communicate()
    if process.returncode:
        raise RuntimeError(f"command failed: {command}\n{stdout}\n{stderr}")
    payload = json.loads(output.read_text(encoding="utf-8"))
    payload["process_wall_seconds"] = time.perf_counter() - started
    payload["sampled_peak_rss_bytes"] = peak
    return payload


def stats(rows, field):
    values = [row[field] for row in rows]
    return {
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
        "standard_deviation": statistics.stdev(values) if len(values) > 1 else 0.0,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default=os.sys.executable)
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    if args.repeats < 1:
        raise ValueError("repeats must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = ROOT / "benchmark_data"
    env = os.environ.copy()
    env.update({"OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"})

    with tempfile.TemporaryDirectory(prefix="sct_cpu_bench_") as temporary:
        temporary = Path(temporary)
        commands = {
            "python": lambda out: [args.python, str(ROOT / "benchmarks" / "python_runner.py"), str(data), str(out)],
            "r": lambda out: [args.rscript, str(ROOT / "benchmarks" / "r_runner.R"), str(data), str(out)],
        }
        for name in ("python", "r"):
            out = temporary / f"warmup_{name}.json"
            monitored(commands[name](out), out, env)
        runs = {"python": [], "r": []}
        for index in range(args.repeats):
            order = ("python", "r") if index % 2 == 0 else ("r", "python")
            for name in order:
                out = temporary / f"run_{index}_{name}.json"
                runs[name].append(monitored(commands[name](out), out, env))

        oracle = temporary / "oracle"
        subprocess.run([
            args.rscript,
            str(ROOT / "benchmarks" / "reference.R"),
            str(data),
            str(oracle),
        ], check=True, env=env)
        quality_path = args.output_dir / "quality.json"
        subprocess.run([
            args.python,
            str(ROOT / "benchmarks" / "compare.py"),
            str(data),
            str(oracle),
            str(quality_path),
        ], check=True, env=env)

    fields = ("compute_wall_seconds", "compute_cpu_seconds",
              "process_wall_seconds", "sampled_peak_rss_bytes")
    summary = {
        "method": "one warm-up, fresh subprocesses, alternating implementation order",
        "repeats": args.repeats,
        "thread_environment": {name: env[name] for name in
                               ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")},
        "platform": platform.platform(),
        "physical_cores": psutil.cpu_count(logical=False),
        "logical_cores": psutil.cpu_count(logical=True),
        "input": {"genes": 800, "cells": 240, "nnz": 22131},
        "summary": {
            name: {field: stats(rows, field) for field in fields}
            for name, rows in runs.items()
        },
        "runs": runs,
        "quality_report": "results/quality.json",
        "interpretation_limit": (
            "Small fixed Windows benchmark; sampled RSS includes interpreter/runtime and may miss sub-10ms peaks."
        ),
    }
    output = args.output_dir / "benchmark.json"
    output.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
