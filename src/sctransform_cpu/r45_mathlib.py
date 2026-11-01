"""Positive-argument R 4.5 digamma/trigamma research compatibility.

This is a narrow Python translation of the ``n in {0, 1}, kode=1, m=1``
path through R 4.5's ``src/nmath/polygamma.c`` (GPL-2-or-later).  It exists to
measure whether R Mathlib primitive differences explain glmGamPoi alignment;
it is not a general psigamma implementation.
"""
from __future__ import annotations

import math
import sys
from decimal import Decimal, localcontext
from functools import lru_cache

import numpy as np

@lru_cache(maxsize=4096)
def exp(value: float) -> float:
    """High-precision exponential for the pinned compatibility candidate.

    MSVC's scalar ``exp`` differs by one ULP from the Rtools45 build at a
    traced PORT iterate that changes the returned theta.  Decimal evaluation
    reproduces that iterate and all four pinned end-to-end regressions.  It is
    intentionally used only by the strict candidate; it is not a claim that
    every R platform or every R-language ``exp`` call is correctly rounded.
    """
    value = float(value)
    if math.isnan(value):
        return math.nan
    if value == math.inf:
        return math.inf
    if value == -math.inf:
        return 0.0
    with localcontext() as context:
        context.prec = 60
        return float(Decimal.from_float(value).exp())


_BERNOULLI = (
    1.00000000000000000e00,
    -5.00000000000000000e-01,
    1.66666666666666667e-01,
    -3.33333333333333333e-02,
    2.38095238095238095e-02,
    -3.33333333333333333e-02,
    7.57575757575757576e-02,
    -2.53113553113553114e-01,
    1.16666666666666667e00,
    -7.09215686274509804e00,
    5.49711779448621554e01,
    -5.29124242424242424e02,
    6.19212318840579710e03,
    -8.65802531135531136e04,
    1.42551716666666667e06,
    -2.72982310678160920e07,
    6.01580873900642368e08,
    -1.51163157670921569e10,
    4.29614643061166667e11,
    -1.37116552050883328e13,
    4.88332318973593167e14,
    -1.92965793419400681e16,
)

_GAMMA_COEFFICIENTS = (
    +.8571195590989331421920062399942e-2, +.4415381324841006757191315771652e-2,
    +.5685043681599363378632664588789e-1, -.4219835396418560501012500186624e-2,
    +.1326808181212460220584006796352e-2, -.1893024529798880432523947023886e-3,
    +.3606925327441245256578082217225e-4, -.6056761904460864218485548290365e-5,
    +.1055829546302283344731823509093e-5, -.1811967365542384048291855891166e-6,
    +.3117724964715322277790254593169e-7, -.5354219639019687140874081024347e-8,
    +.9193275519859588946887786825940e-9, -.1577941280288339761767423273953e-9,
    +.2707980622934954543266540433089e-10, -.4646818653825730144081661058933e-11,
    +.7973350192007419656460767175359e-12, -.1368078209830916025799499172309e-12,
    +.2347319486563800657233471771688e-13, -.4027432614949066932766570534699e-14,
    +.6910051747372100912138336975257e-15, -.1185584500221992907052387126192e-15,
)

_LGAMMA_CORRECTION = (
    +.1666389480451863247205729650822e0,
    -.1384948176067563840732986059135e-4,
    +.9810825646924729426157171547487e-8,
    -.1809129475572494194263306266719e-10,
    +.6221098041892605227126015543416e-13,
)


def _chebyshev_eval(x: float, coefficients) -> float:
    twox = x * 2.0
    b0 = b1 = b2 = 0.0
    for coefficient in reversed(coefficients):
        b2, b1 = b1, b0
        b0 = twox * b1 - b2 + coefficient
    return (b0 - b2) * 0.5


def _gamma_positive(x: float) -> float:
    n = int(x)
    y = x - n
    n -= 1
    value = _chebyshev_eval(y * 2.0 - 1.0, _GAMMA_COEFFICIENTS) + 0.9375
    if n == 0:
        return value
    if n < 0:
        for i in range(-n):
            value /= x + i
        return value
    for i in range(1, n + 1):
        value *= y + i
    return value


def _lgammacor(x: float) -> float:
    if x < 94_906_265.62425156:
        temporary = 10.0 / x
        return _chebyshev_eval(temporary * temporary * 2.0 - 1.0,
                               _LGAMMA_CORRECTION) / x
    return 1.0 / (x * 12.0)


def _lgamma_positive_single(x: float) -> float:
    x = float(x)
    if not math.isfinite(x) or x <= 0.0:
        raise ValueError("the compatibility subset requires finite x>0")
    if x < 1e-306:
        return -math.log(x)
    if x <= 10.0:
        return math.log(abs(_gamma_positive(x)))
    if x > 1e17:
        return x * (math.log(x) - 1.0)
    log_sqrt_two_pi = 0.91893853320467274178032973640562
    if x > 4_934_720.0:
        return log_sqrt_two_pi + (x - 0.5) * math.log(x) - x
    return log_sqrt_two_pi + (x - 0.5) * math.log(x) - x + _lgammacor(x)


def _dpsifn_positive_single(x: float, n: int) -> float:
    """Return R's scaled dpsifn result for positive x and n=0 or n=1."""
    x = float(x)
    if not math.isfinite(x) or x <= 0.0 or n not in (0, 1):
        raise ValueError("the compatibility subset requires finite x>0 and n in {0,1}")
    xln = math.log(x)
    large = 1.0 / (2.0 * sys.float_info.epsilon)
    if n == 0 and x * xln > large:
        return -xln
    if n >= 1 and x > n * large:
        return math.exp(-n * xln) / n

    # R's i1mach/d1mach values on IEEE-754 binary64.
    nx = 1021
    r1m5 = 0.301029995663981195213738894724493027
    r1m4 = sys.float_info.epsilon * 0.5
    wdtol = max(r1m4, 0.5e-18)
    elim = 2.302 * (nx * r1m5 - 3.0)
    rln = min(r1m5 * 53.0, 18.06)
    fln = max(rln, 3.0) - 3.0
    yint = 3.50 + 0.40 * fln
    slope = 0.21 + fln * (0.0006038 * fln + 0.008677)

    nn = n
    fn = float(nn)
    t = (fn + 1.0) * xln
    if abs(t) <= elim:
        if x < wdtol:
            return x ** (-n - 1)
        xm = yint + slope * fn
        xmin = float(int(xm) + 1)
        use_series = False
        if n != 0:
            xm = -2.302 * rln - min(0.0, xln)
            arg = min(0.0, xm / n)
            eps = math.exp(arg)
            xm = -arg if abs(arg) < 1.0e-3 else 1.0 - eps
            fln = x * xm / eps
            xm = xmin - x
            use_series = xm > 7.0 and fln < 15.0
        if use_series:
            terms = int(fln) + 1
            npower = n + 1
            value = math.exp(-(n + 1.0) * xln)
            total = value
            denominator = x
            for _ in range(1, terms + 1):
                denominator += 1.0
                total += math.pow(denominator, float(-npower))
            return total

        xdmy = x
        xdmln = xln
        xinc = 0.0
        if x < xmin:
            integer_x = int(x)
            xinc = xmin - integer_x
            xdmy = x + xinc
            xdmln = math.log(xdmy)

        t = fn * xdmln
        t1 = xdmln + xdmln
        t2 = t + xdmln
        if max(abs(t), abs(t1), abs(t2)) <= elim:
            tss = math.exp(-t)
            tt = 0.5 / xdmy
            t1 = tt
            tst = wdtol * tt
            if nn != 0:
                t1 = tt + 1.0 / fn
            rxsq = 1.0 / (xdmy * xdmy)
            ta = 0.5 * rxsq
            t = (fn + 1.0) * ta
            total = t * _BERNOULLI[2]
            if abs(total) >= tst:
                tk = 2.0
                for k in range(4, 23):
                    t = (
                        t
                        * ((tk + fn + 1.0) / (tk + 1.0))
                        * ((tk + fn) / (tk + 2.0))
                        * rxsq
                    )
                    term = t * _BERNOULLI[k - 1]
                    if abs(term) < tst:
                        break
                    total += term
                    tk += 2.0
            total = (total + t1) * tss
            if xinc != 0.0:
                recurrence_count = int(xinc)
                if nn == 0:
                    for i in range(1, recurrence_count + 1):
                        total += 1.0 / (x + (recurrence_count - i))
                else:
                    power = nn + 1
                    xm = xinc - 1.0
                    fx = x + xm
                    for _ in range(1, recurrence_count + 1):
                        total += math.pow(fx, float(-power))
                        xm -= 1.0
                        fx = x + xm
            return total - xdmln if nn == 0 else total

    raise OverflowError("dpsifn compatibility subset overflowed")


def _apply(values, n: int):
    array = np.asarray(values, dtype=np.float64)
    result = np.empty_like(array)
    source = array.ravel()
    target = result.ravel()
    for i, value in enumerate(source):
        target[i] = _dpsifn_positive_single(float(value), n)
    return float(result) if result.ndim == 0 else result


def digamma(values):
    """R 4.5 digamma for finite positive binary64 inputs."""
    return -_apply(values, 0)


def trigamma(values):
    """R 4.5 trigamma for finite positive binary64 inputs."""
    return _apply(values, 1)


def lgamma(values):
    """R 4.5 lgammafn for finite positive binary64 inputs."""
    array = np.asarray(values, dtype=np.float64)
    result = np.empty_like(array)
    for i, value in enumerate(array.ravel()):
        result.ravel()[i] = _lgamma_positive_single(float(value))
    return float(result) if result.ndim == 0 else result
