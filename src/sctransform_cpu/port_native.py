"""Pinned R 4.5.2 PORT compatibility backend.

The bundled native library links R's unmodified ``portsrc.f`` reverse-
communication routines to a deliberately tiny C adapter.  Python still owns
the statistical objective, gradient, and Hessian; the native component only
reproduces the optimizer state machine and R/MinGW extended-precision start
statistics needed for strict oracle parity.

R is GPL-2-or-later.  The adapter and distributions containing the linked
library are GPL-3.0-or-later.  This module fails closed when the pinned native
library is unavailable; it never silently substitutes another optimizer.
"""
from __future__ import annotations

from dataclasses import dataclass
import ctypes
import math
import os
from pathlib import Path
from typing import Callable

import numpy as np


_MESSAGES = {
    3: "X-convergence (3)",
    4: "relative convergence (4)",
    5: "both X-convergence and relative convergence (5)",
    6: "absolute function convergence (6)",
    7: "singular convergence (7)",
    8: "false convergence (8)",
    9: "function evaluation limit reached without convergence (9)",
    10: "iteration limit reached without convergence (10)",
    14: "storage only has been allocated (14)",
    15: "LIV too small (15)",
    16: "LV too small (16)",
    63: "objective cannot be computed at initial parameter (63)",
    65: "gradient cannot be computed at initial parameter (65)",
    300: "initial parameter violates constraints (300)",
}


@dataclass(frozen=True)
class NativePortResult:
    x: float
    objective: float
    gradient: float
    hessian: float | None
    status: int
    convergence: int
    message: str
    iterations: int
    function_evaluations: int
    gradient_evaluations: int


_LIBRARY: ctypes.CDLL | None = None
_LIBRARY_PATH: Path | None = None


def _candidate_library_paths() -> list[Path]:
    override = os.environ.get("SCT_PORT_R45_DLL")
    package = Path(__file__).resolve().parent
    candidates = []
    if override:
        candidates.append(Path(override))
    candidates.extend([
        package / "native" / "sct_port_r45.dll",
        package.parent / "build" / "sct_port_r45.dll",
    ])
    return candidates


def _load_library() -> ctypes.CDLL:
    global _LIBRARY, _LIBRARY_PATH
    if _LIBRARY is not None:
        return _LIBRARY
    path = next((candidate for candidate in _candidate_library_paths()
                 if candidate.is_file()), None)
    if path is None:
        searched = ", ".join(str(item) for item in _candidate_library_paths())
        raise RuntimeError(
            "strict R 4.5.2 PORT backend is unavailable; searched: " + searched
        )
    library = ctypes.CDLL(str(path))
    double_pointer = ctypes.POINTER(ctypes.c_double)
    int_pointer = ctypes.POINTER(ctypes.c_int32)
    library.sct_port_state_size.argtypes = []
    library.sct_port_state_size.restype = ctypes.c_int
    library.sct_port_initialize.argtypes = [
        int_pointer, ctypes.c_int, double_pointer, ctypes.c_int, ctypes.c_int,
    ]
    library.sct_port_initialize.restype = None
    library.sct_port_iterate_bfgs.argtypes = [
        double_pointer, double_pointer, double_pointer, double_pointer,
        int_pointer, ctypes.c_int, ctypes.c_int, double_pointer, double_pointer,
    ]
    library.sct_port_iterate_bfgs.restype = None
    library.sct_port_iterate_hessian.argtypes = [
        double_pointer, double_pointer, double_pointer, double_pointer,
        double_pointer, int_pointer, ctypes.c_int, ctypes.c_int,
        double_pointer, double_pointer,
    ]
    library.sct_port_iterate_hessian.restype = None
    callback_arguments = [
        double_pointer, double_pointer, ctypes.c_int, ctypes.c_double,
        ctypes.c_int, double_pointer, double_pointer, ctypes.c_int,
    ]
    library.sct_glmgampoi_score.argtypes = callback_arguments
    library.sct_glmgampoi_score.restype = ctypes.c_double
    library.sct_glmgampoi_score_derivative.argtypes = callback_arguments
    library.sct_glmgampoi_score_derivative.restype = ctypes.c_double
    _LIBRARY, _LIBRARY_PATH = library, path.resolve()
    return library


def native_library_path() -> Path:
    """Return the loaded strict backend path, loading it if necessary."""
    _load_library()
    assert _LIBRARY_PATH is not None
    return _LIBRARY_PATH


def _native_callback_arguments(y, mu, alpha, cox_reid, table):
    y_array = np.ascontiguousarray(y, dtype=np.float64)
    mu_array = np.ascontiguousarray(mu, dtype=np.float64)
    if y_array.ndim != 1 or mu_array.shape != y_array.shape:
        raise ValueError("expected equal one-dimensional y and mu")
    if not math.isfinite(alpha) or alpha <= 0.0:
        raise ValueError("alpha must be finite and positive")
    if table is None:
        keys = frequencies = np.empty(0, dtype=np.float64)
    else:
        keys = np.ascontiguousarray(table[0], dtype=np.float64)
        frequencies = np.ascontiguousarray(table[1], dtype=np.float64)
        if keys.ndim != 1 or frequencies.shape != keys.shape:
            raise ValueError("invalid count table")
    pointer = ctypes.POINTER(ctypes.c_double)
    null = pointer()
    arguments = (
        y_array.ctypes.data_as(pointer),
        mu_array.ctypes.data_as(pointer),
        len(y_array),
        float(alpha),
        int(bool(cox_reid)),
        keys.ctypes.data_as(pointer) if len(keys) else null,
        frequencies.ctypes.data_as(pointer) if len(frequencies) else null,
        len(keys),
    )
    # Keep the arrays alive until the ctypes call has returned.
    return arguments, (y_array, mu_array, keys, frequencies)


def glmgampoi_score_native(y, mu, alpha, *, cox_reid=True, table=None) -> float:
    """Evaluate the pinned Rtools scalar glmGamPoi score callback."""
    library = _load_library()
    arguments, keepalive = _native_callback_arguments(
        y, mu, alpha, cox_reid, table)
    value = float(library.sct_glmgampoi_score(*arguments))
    del keepalive
    return value


def glmgampoi_score_derivative_native(
    y, mu, alpha, *, cox_reid=True, table=None
) -> float:
    """Evaluate the pinned Rtools scalar glmGamPoi score derivative."""
    library = _load_library()
    arguments, keepalive = _native_callback_arguments(
        y, mu, alpha, cox_reid, table)
    value = float(library.sct_glmgampoi_score_derivative(*arguments))
    del keepalive
    return value


def minimize_port_native_1d(
    objective: Callable[[float], float],
    gradient: Callable[[float], float],
    start: float,
    lower: float,
    upper: float,
    *,
    hessian: Callable[[float], float] | None = None,
    max_iter: int = 200,
) -> NativePortResult:
    """Run the pinned one-dimensional bounded PORT state machine.

    ``hessian=None`` selects R's gradient-only ``DRMNGB`` path; providing a
    Hessian selects ``DRMNHB``, exactly as ``stats::nlminb`` does.
    """
    if not (math.isfinite(start) and lower <= start <= upper and lower < upper):
        raise ValueError("start and bounds must be finite and ordered")
    if max_iter <= 0:
        raise ValueError("max_iter must be positive")

    library = _load_library()
    liv = max(int(library.sct_port_state_size()), 128)
    lv = 256
    iv = np.zeros(liv, dtype=np.int32)
    v = np.zeros(lv, dtype=np.float64)
    bounds = np.asarray([lower, upper], dtype=np.float64)
    scale = np.asarray([1.0], dtype=np.float64)
    x = np.asarray([start], dtype=np.float64)
    g = np.zeros(1, dtype=np.float64)
    h = np.zeros(1, dtype=np.float64)
    fx = ctypes.c_double(0.0)

    int_pointer = ctypes.POINTER(ctypes.c_int32)
    double_pointer = ctypes.POINTER(ctypes.c_double)
    iv_pointer = iv.ctypes.data_as(int_pointer)
    v_pointer = v.ctypes.data_as(double_pointer)
    bounds_pointer = bounds.ctypes.data_as(double_pointer)
    scale_pointer = scale.ctypes.data_as(double_pointer)
    x_pointer = x.ctypes.data_as(double_pointer)
    g_pointer = g.ctypes.data_as(double_pointer)
    h_pointer = h.ctypes.data_as(double_pointer)
    library.sct_port_initialize(iv_pointer, liv, v_pointer, lv, max_iter)

    # Every iteration advances one reverse-communication event.  The guard is
    # intentionally much larger than PORT's evaluation limits but finite, so
    # a malformed native state cannot hang the Python process indefinitely.
    for _ in range(20_000):
        if hessian is None:
            library.sct_port_iterate_bfgs(
                bounds_pointer, scale_pointer, ctypes.byref(fx), g_pointer,
                iv_pointer, liv, lv, v_pointer, x_pointer,
            )
        else:
            library.sct_port_iterate_hessian(
                bounds_pointer, scale_pointer, ctypes.byref(fx), g_pointer,
                h_pointer, iv_pointer, liv, lv, v_pointer, x_pointer,
            )
        status = int(iv[0])
        if status == 1:
            fx.value = float(objective(float(x[0])))
            if not math.isfinite(fx.value):
                iv[0] = 63
        elif status == 2:
            g[0] = float(gradient(float(x[0])))
            if hessian is not None:
                h[0] = float(hessian(float(x[0])))
            if not math.isfinite(float(g[0])) or (
                hessian is not None and not math.isfinite(float(h[0]))
            ):
                iv[0] = 65
        else:
            break
    else:
        raise RuntimeError("PORT reverse-communication guard exhausted")

    status = int(iv[0])
    convergence = 0 if status in range(3, 7) else 1
    return NativePortResult(
        x=float(x[0]),
        objective=float(v[9]),
        gradient=float(g[0]),
        hessian=float(h[0]) if hessian is not None else None,
        status=status,
        convergence=convergence,
        message=_MESSAGES.get(status, f"PORT status ({status})"),
        iterations=int(iv[30]),
        function_evaluations=int(iv[5]),
        gradient_evaluations=int(iv[29]),
    )
