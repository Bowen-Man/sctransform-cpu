from __future__ import annotations

import os
from pathlib import Path
import subprocess

import pytest

from benchmarks.bit_identity import compare_oracle


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
