"""Pinned R 4.5.2 PORT compatibility backend.

The bundled native library links R's unmodified ``portsrc.f`` reverse-
communication routines to narrow C/C++ adapters. Python owns orchestration;
the native boundary reproduces the optimizer state machine and the exact
R/Rcpp, glmGamPoi and matrixStats arithmetic needed for strict oracle parity.

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
import threading
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
_FPU_THREAD = threading.local()


def _ensure_r_fpu(library: ctypes.CDLL) -> None:
    """Set R's x87 precision once for each Python thread using the backend."""
    if not getattr(_FPU_THREAD, "configured", False):
        library.sct_use_r_fpu()
        _FPU_THREAD.configured = True


def _candidate_library_paths() -> list[Path]:
    override = (os.environ.get("SCT_PORT_NATIVE_LIBRARY")
                or os.environ.get("SCT_PORT_R45_DLL"))
    package = Path(__file__).resolve().parent
    native_name = ("sct_port_r45.dll" if os.name == "nt"
                   else "libsct_port_r45.so")
    candidates = []
    if override:
        candidates.append(Path(override))
    candidates.extend([
        package / "native" / native_name,
        package.parents[1] / "build" / "native-linux" / native_name,
    ])
    return candidates


def _load_library() -> ctypes.CDLL:
    global _LIBRARY, _LIBRARY_PATH
    if _LIBRARY is not None:
        _ensure_r_fpu(_LIBRARY)
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
    library.sct_glmgampoi_loglikelihood.argtypes = callback_arguments
    library.sct_glmgampoi_loglikelihood.restype = ctypes.c_double
    library.sct_glmgampoi_score_log.argtypes = callback_arguments
    library.sct_glmgampoi_score_log.restype = ctypes.c_double
    library.sct_glmgampoi_score_derivative_log.argtypes = callback_arguments
    library.sct_glmgampoi_score_derivative_log.restype = ctypes.c_double
    library.sct_r45_exp.argtypes = [ctypes.c_double]
    library.sct_r45_exp.restype = ctypes.c_double
    library.sct_use_r_fpu.argtypes = []
    library.sct_use_r_fpu.restype = None
    library.sct_row_gmean.argtypes = [
        double_pointer, int_pointer, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        double_pointer, double_pointer,
    ]
    library.sct_row_gmean.restype = None
    library.sct_fit_beta_one_group.argtypes = [
        double_pointer, double_pointer, ctypes.c_int,
        ctypes.c_double, ctypes.c_double,
    ]
    library.sct_fit_beta_one_group.restype = ctypes.c_double
    library.sct_make_mu.argtypes = [
        ctypes.c_double, double_pointer, ctypes.c_int, double_pointer,
    ]
    library.sct_make_mu.restype = None
    library.sct_cpp_exp.argtypes = [ctypes.c_double]
    library.sct_cpp_exp.restype = ctypes.c_double
    library.sct_make_count_table.argtypes = [
        double_pointer, ctypes.c_int, ctypes.c_int,
        double_pointer, double_pointer,
    ]
    library.sct_make_count_table.restype = ctypes.c_int
    library.sct_count_mean_variance.argtypes = [
        double_pointer, ctypes.c_int, double_pointer,
    ]
    library.sct_count_mean_variance.restype = None
    if os.name != "nt":
        library.sct_log10_array.argtypes = [
            double_pointer, ctypes.c_int, double_pointer,
        ]
        library.sct_log10_array.restype = None
        library.sct_reciprocal_array.argtypes = [
            double_pointer, ctypes.c_int, double_pointer,
        ]
        library.sct_reciprocal_array.restype = None
        library.sct_regularization_dispersion.argtypes = [
            double_pointer, double_pointer, ctypes.c_int, double_pointer,
        ]
        library.sct_regularization_dispersion.restype = None
        library.sct_regularization_theta.argtypes = [
            double_pointer, double_pointer, ctypes.c_int, double_pointer,
        ]
        library.sct_regularization_theta.restype = None
    library.sct_ksmooth_normal.argtypes = [
        double_pointer, double_pointer, ctypes.c_int,
        double_pointer, ctypes.c_int, ctypes.c_double, double_pointer,
    ]
    library.sct_ksmooth_normal.restype = None
    library.sct_bw_sj.argtypes = [double_pointer, ctypes.c_int, ctypes.c_int]
    library.sct_bw_sj.restype = ctypes.c_double
    library.sct_pearson_residuals.argtypes = [
        double_pointer, ctypes.c_int, ctypes.c_int, double_pointer,
        double_pointer, double_pointer, double_pointer,
        ctypes.c_double, ctypes.c_double, double_pointer,
        double_pointer, double_pointer, double_pointer,
    ]
    library.sct_pearson_residuals.restype = None
    library.sct_residual_statistics.argtypes = [
        double_pointer, ctypes.c_int, ctypes.c_int,
        double_pointer, double_pointer,
    ]
    library.sct_residual_statistics.restype = None
    # R 4.5.2 initializes the x87 precision field to extended precision
    # (control word 0x037f), whereas CPython on Windows uses 0x027f.  The
    # difference is observable in Rcpp's exp(x)-1 and glmGamPoi callbacks.
    # Establish the pinned R floating-point environment on this thread before
    # any strict compatibility arithmetic runs.
    _ensure_r_fpu(library)
    _LIBRARY, _LIBRARY_PATH = library, path.resolve()
    return library


def native_library_path() -> Path:
    """Return the loaded strict backend path, loading it if necessary."""
    _load_library()
    assert _LIBRARY_PATH is not None
    return _LIBRARY_PATH


def r45_exp_native(value: float) -> float:
    """Return the pinned Rtools45 scalar exponential."""
    return float(_load_library().sct_r45_exp(float(value)))


def count_table_native(values, stop_if_larger: int):
    """Build glmGamPoi's count table in native unordered-map order."""
    array = np.ascontiguousarray(values, dtype=np.float64)
    if array.ndim != 1 or stop_if_larger < 1:
        raise ValueError("count table input must be a vector with a positive limit")
    keys = np.empty(stop_if_larger + 1, dtype=np.float64)
    frequencies = np.empty(stop_if_larger + 1, dtype=np.float64)
    pointer = ctypes.POINTER(ctypes.c_double)
    size = int(_load_library().sct_make_count_table(
        array.ctypes.data_as(pointer),
        len(array),
        stop_if_larger,
        keys.ctypes.data_as(pointer),
        frequencies.ctypes.data_as(pointer),
    ))
    if size == 0:
        return None
    return keys[:size].copy(), frequencies[:size].copy()


def platform_log10_native(values) -> np.ndarray:
    """Evaluate scalar log10 with the pinned platform C runtime."""
    array = np.ascontiguousarray(values, dtype=np.float64)
    if array.ndim != 1 or np.any(array <= 0.0):
        raise ValueError("log10 input must be a positive one-dimensional vector")
    if os.name == "nt":
        return np.log10(array)
    result = np.empty_like(array)
    pointer = ctypes.POINTER(ctypes.c_double)
    _load_library().sct_log10_array(
        array.ctypes.data_as(pointer), len(array),
        result.ctypes.data_as(pointer))
    return result


def reciprocal_native(values) -> np.ndarray:
    """Evaluate the scalar R/C reciprocal loop without NumPy reassociation."""
    array = np.ascontiguousarray(values, dtype=np.float64)
    if array.ndim != 1:
        raise ValueError("reciprocal input must be a one-dimensional vector")
    if os.name == "nt":
        return np.divide(
            1.0, array, out=np.full_like(array, np.inf), where=array > 0)
    result = np.empty_like(array)
    pointer = ctypes.POINTER(ctypes.c_double)
    _load_library().sct_reciprocal_array(
        array.ctypes.data_as(pointer), len(array),
        result.ctypes.data_as(pointer))
    return result


def regularization_dispersion_native(log_gmean, theta) -> np.ndarray:
    """Match R's ``log10(1 + 10^x / theta)`` scalar evaluation."""
    x = np.ascontiguousarray(log_gmean, dtype=np.float64)
    theta_array = np.ascontiguousarray(theta, dtype=np.float64)
    if x.ndim != 1 or theta_array.shape != x.shape:
        raise ValueError("regularization inputs must be aligned vectors")
    if os.name == "nt":
        return np.log10(1.0 + np.power(10.0, x) / theta_array)
    result = np.empty_like(x)
    pointer = ctypes.POINTER(ctypes.c_double)
    _load_library().sct_regularization_dispersion(
        x.ctypes.data_as(pointer), theta_array.ctypes.data_as(pointer), len(x),
        result.ctypes.data_as(pointer))
    return result


def regularization_theta_native(log_gmean, dispersion) -> np.ndarray:
    """Match R's ``10^x / (10^dispersion - 1)`` scalar evaluation."""
    x = np.ascontiguousarray(log_gmean, dtype=np.float64)
    dispersion_array = np.ascontiguousarray(dispersion, dtype=np.float64)
    if x.ndim != 1 or dispersion_array.shape != x.shape:
        raise ValueError("regularization inputs must be aligned vectors")
    if os.name == "nt":
        return np.power(10.0, x) / (np.power(10.0, dispersion_array) - 1.0)
    result = np.empty_like(x)
    pointer = ctypes.POINTER(ctypes.c_double)
    _load_library().sct_regularization_theta(
        x.ctypes.data_as(pointer), dispersion_array.ctypes.data_as(pointer),
        len(x), result.ctypes.data_as(pointer))
    return result


def row_gmean_native(counts) -> tuple[np.ndarray, np.ndarray]:
    """Return pinned row geometric means and their R-style base-10 logs."""
    from scipy import sparse

    matrix = sparse.csc_matrix(counts, dtype=np.float64, copy=True)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    matrix.sort_indices()
    values = np.ascontiguousarray(matrix.data, dtype=np.float64)
    row_index = np.ascontiguousarray(matrix.indices, dtype=np.int32)
    rows, cols = matrix.shape
    if cols <= 0:
        raise ValueError("matrix must contain cells")
    gmean = np.empty(rows, dtype=np.float64)
    log_gmean = np.empty(rows, dtype=np.float64)
    double_pointer = ctypes.POINTER(ctypes.c_double)
    int_pointer = ctypes.POINTER(ctypes.c_int32)
    _load_library().sct_row_gmean(
        values.ctypes.data_as(double_pointer),
        row_index.ctypes.data_as(int_pointer),
        len(values),
        rows,
        cols,
        gmean.ctypes.data_as(double_pointer),
        log_gmean.ctypes.data_as(double_pointer),
    )
    if os.name != "nt":
        log_gmean = platform_log10_native(gmean)
    return gmean, log_gmean


def fit_beta_one_group_native(y, offset, alpha: float, initial: float) -> float:
    """Run glmGamPoi 1.22's one-group Newton loop in the pinned ABI."""
    counts = np.ascontiguousarray(y, dtype=np.float64)
    offsets = np.ascontiguousarray(offset, dtype=np.float64)
    if counts.ndim != 1 or offsets.shape != counts.shape or len(counts) == 0:
        raise ValueError("expected aligned nonempty vectors")
    pointer = ctypes.POINTER(ctypes.c_double)
    return float(_load_library().sct_fit_beta_one_group(
        counts.ctypes.data_as(pointer), offsets.ctypes.data_as(pointer),
        len(counts), float(alpha), float(initial)))


def make_mu_native(beta: float, offset) -> np.ndarray:
    """Evaluate glmGamPoi's exp(beta + offset) expression exactly."""
    offsets = np.ascontiguousarray(offset, dtype=np.float64)
    if offsets.ndim != 1:
        raise ValueError("expected a one-dimensional offset")
    result = np.empty_like(offsets)
    pointer = ctypes.POINTER(ctypes.c_double)
    _load_library().sct_make_mu(
        float(beta), offsets.ctypes.data_as(pointer), len(offsets),
        result.ctypes.data_as(pointer))
    return result


def cpp_exp_native(value: float) -> float:
    """Evaluate the pinned Rcpp scalar exponential."""
    return float(_load_library().sct_cpp_exp(float(value)))


def count_mean_variance_native(values) -> tuple[float, float]:
    """Reproduce R mean/var reductions under its x87 long-double ABI."""
    array = np.ascontiguousarray(values, dtype=np.float64)
    if array.ndim != 1 or len(array) < 2:
        raise ValueError("sample variance requires at least two values")
    result = np.empty(4, dtype=np.float64)
    pointer = ctypes.POINTER(ctypes.c_double)
    _load_library().sct_count_mean_variance(
        array.ctypes.data_as(pointer), len(array),
        result.ctypes.data_as(pointer))
    return float(result[0]), float(result[1])


def count_dispersion_start_native(values) -> tuple[float, float]:
    """Return glmGamPoi's R-evaluated moment start and its natural log."""
    array = np.ascontiguousarray(values, dtype=np.float64)
    if array.ndim != 1 or len(array) < 2:
        raise ValueError("dispersion start requires at least two values")
    result = np.empty(4, dtype=np.float64)
    pointer = ctypes.POINTER(ctypes.c_double)
    _load_library().sct_count_mean_variance(
        array.ctypes.data_as(pointer), len(array),
        result.ctypes.data_as(pointer))
    start = float(result[2])
    if not math.isfinite(start) or start <= 0.0:
        return 0.5, math.log(0.5)
    return start, float(result[3])


def ksmooth_normal_native(x, y, query, bandwidth: float) -> np.ndarray:
    """Run R 4.5.2's normal-kernel BDRksmooth loop exactly."""
    x_array = np.asarray(x, dtype=np.float64)
    y_array = np.asarray(y, dtype=np.float64)
    query_array = np.asarray(query, dtype=np.float64)
    if (x_array.ndim != 1 or y_array.shape != x_array.shape
            or query_array.ndim != 1 or len(x_array) == 0
            or len(query_array) == 0 or not math.isfinite(bandwidth)
            or bandwidth <= 0.0):
        raise ValueError("invalid ksmooth input")
    order = np.argsort(x_array, kind="stable")
    query_order = np.argsort(query_array, kind="stable")
    sorted_x = np.ascontiguousarray(x_array[order])
    sorted_y = np.ascontiguousarray(y_array[order])
    sorted_query = np.ascontiguousarray(query_array[query_order])
    sorted_result = np.empty_like(sorted_query)
    pointer = ctypes.POINTER(ctypes.c_double)
    _load_library().sct_ksmooth_normal(
        sorted_x.ctypes.data_as(pointer), sorted_y.ctypes.data_as(pointer),
        len(sorted_x), sorted_query.ctypes.data_as(pointer), len(sorted_query),
        float(bandwidth), sorted_result.ctypes.data_as(pointer))
    result = np.empty_like(sorted_result)
    result[query_order] = sorted_result
    return result


def bw_sj_native(values, bins: int = 1000) -> float:
    """Run the pinned R bw.SJ STE path, including its zeroin routine."""
    array = np.sort(np.asarray(values, dtype=np.float64), kind="stable")
    if (array.ndim != 1 or len(array) < 2 or not np.all(np.isfinite(array))
            or bins < 1):
        raise ValueError("bw.SJ requires at least two finite values")
    contiguous = np.ascontiguousarray(array)
    pointer = ctypes.POINTER(ctypes.c_double)
    result = float(_load_library().sct_bw_sj(
        contiguous.ctypes.data_as(pointer), len(contiguous), bins))
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("bw.SJ failed")
    return result


def pearson_residuals_native(
    counts, library, intercept, slope, theta, minimum_variance, clip,
    *, return_components: bool = False,
) -> np.ndarray | tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate the pinned SCT Pearson-residual scalar expressions."""
    matrix = np.ascontiguousarray(counts, dtype=np.float64)
    library_array = np.asarray(library, dtype=np.float64)
    intercept_array = np.ascontiguousarray(intercept, dtype=np.float64)
    slope_array = np.ascontiguousarray(slope, dtype=np.float64)
    theta_array = np.ascontiguousarray(theta, dtype=np.float64)
    if (matrix.ndim != 2 or library_array.shape != (matrix.shape[1],)
            or intercept_array.shape != (matrix.shape[0],)
            or slope_array.shape != intercept_array.shape
            or theta_array.shape != intercept_array.shape):
        raise ValueError("invalid residual tile shapes")
    log_umi = platform_log10_native(library_array)
    result = np.empty_like(matrix)
    pointer = ctypes.POINTER(ctypes.c_double)
    mu = np.empty_like(matrix) if return_components else None
    variance = np.empty_like(matrix) if return_components else None
    denominator = np.empty_like(matrix) if return_components else None
    null = pointer()
    _load_library().sct_pearson_residuals(
        matrix.ctypes.data_as(pointer), matrix.shape[0], matrix.shape[1],
        log_umi.ctypes.data_as(pointer),
        intercept_array.ctypes.data_as(pointer),
        slope_array.ctypes.data_as(pointer),
        theta_array.ctypes.data_as(pointer), float(minimum_variance),
        float(clip), result.ctypes.data_as(pointer),
        mu.ctypes.data_as(pointer) if mu is not None else null,
        variance.ctypes.data_as(pointer) if variance is not None else null,
        denominator.ctypes.data_as(pointer) if denominator is not None else null)
    return (result, mu, variance, denominator) if return_components else result


def residual_statistics_native(residual) -> tuple[np.ndarray, np.ndarray]:
    """Match base::rowMeans and matrixStats::rowVars for a dense SCT tile."""
    matrix = np.ascontiguousarray(residual, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[1] < 2:
        raise ValueError("residual tile must be two-dimensional with >=2 cells")
    means = np.empty(matrix.shape[0], dtype=np.float64)
    variances = np.empty(matrix.shape[0], dtype=np.float64)
    pointer = ctypes.POINTER(ctypes.c_double)
    _load_library().sct_residual_statistics(
        matrix.ctypes.data_as(pointer), matrix.shape[0], matrix.shape[1],
        means.ctypes.data_as(pointer), variances.ctypes.data_as(pointer))
    return means, variances


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


def _glmgampoi_log_callback(function_name, y, mu, log_alpha, cox_reid, table):
    library = _load_library()
    arguments, keepalive = _native_callback_arguments(
        y, mu, 1.0, cox_reid, table)
    mutable = list(arguments)
    mutable[3] = float(log_alpha)
    value = float(getattr(library, function_name)(*mutable))
    del keepalive
    return value


def glmgampoi_loglikelihood_native(
    y, mu, log_alpha, *, cox_reid=True, table=None
) -> float:
    return _glmgampoi_log_callback(
        "sct_glmgampoi_loglikelihood", y, mu, log_alpha, cox_reid, table)


def glmgampoi_score_log_native(
    y, mu, log_alpha, *, cox_reid=True, table=None
) -> float:
    return _glmgampoi_log_callback(
        "sct_glmgampoi_score_log", y, mu, log_alpha, cox_reid, table)


def glmgampoi_score_derivative_log_native(
    y, mu, log_alpha, *, cox_reid=True, table=None
) -> float:
    return _glmgampoi_log_callback(
        "sct_glmgampoi_score_derivative_log", y, mu, log_alpha,
        cox_reid, table)


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
