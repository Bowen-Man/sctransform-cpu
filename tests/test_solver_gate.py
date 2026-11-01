from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from sctransform_cpu.port_native import minimize_port_native_1d, native_library_path
from sctransform_cpu.r_glmgampoi import fit_dispersion_port_native


DATA = Path(__file__).parent / "data"


def test_bundled_native_library_and_quadratic():
    assert native_library_path().name == "sct_port_r45.dll"
    result = minimize_port_native_1d(
        lambda value: (value - 2.0) ** 2,
        lambda value: 2.0 * (value - 2.0),
        0.0,
        -10.0,
        10.0,
        hessian=lambda value: 2.0,
    )
    assert result.convergence == 0
    assert abs(result.x - 2.0) <= 1e-14


def test_stored_r452_solver_gate():
    fixture = np.load(DATA / "solver_stress.npz", allow_pickle=False)
    rows = json.loads((DATA / "solver_stress_r452.json").read_text(encoding="utf-8"))["rows"]
    expected = np.asarray([row["r_alpha"] for row in rows])
    actual = np.asarray([
        fit_dispersion_port_native(y, mu)
        for y, mu in zip(fixture["y"], fixture["mu"])
    ])
    allowed = 1e-8 + 1e-5 * np.abs(expected)
    finite = np.isfinite(actual) & np.isfinite(expected)
    assert np.all(finite)
    assert np.count_nonzero(np.abs(actual - expected) > allowed) == 0
    assert np.array_equal(actual == 0.0, expected == 0.0)
