"""glmGamPoi 1.22 single-intercept dispersion compatibility implementation.

The objective, score and score derivative follow
glmGamPoi ``src/overdispersion.cpp`` (GPL-3 package; portions there are marked
as adapted from DESeq2 under LGPL-3-or-later). The public pipeline uses
``fit_dispersion_port_native``; the other solvers remain private diagnostic
references for auditing the native path.
"""
from __future__ import annotations

import math
import os
import numpy as np
from scipy.special import digamma, gammaln, polygamma


_CR = 0.99
_LOWER = math.log(1e-16)
_UPPER = math.log(1e16)


def _table_if_small(y: np.ndarray):
    """Return count/frequency arrays when distinct values occupy at most n/2.

    glmGamPoi 1.22 builds this table with libstdc++ ``unordered_map`` and then
    iterates it directly.  For collision-free integer keys after its reserve,
    the observed GCC order is reverse first-insertion order.  Keeping that
    order matters because large digamma terms nearly cancel.  Collision-heavy
    tables still require broader conformance tests before this is promoted.
    """
    counts: dict[int, int] = {}
    limit = len(y) // 2
    for value in y:
        key = int(value)
        counts[key] = counts.get(key, 0) + 1
        if len(counts) > limit:
            return None
    items = list(counts.items())[::-1]
    return (np.fromiter((item[0] for item in items), dtype=np.float64),
            np.fromiter((item[1] for item in items), dtype=np.float64))


def _sequential_sum(values) -> float:
    result = 0.0
    for value in values:
        result += float(value)
    return result


def _armadillo_generic_dot_ones(values) -> float:
    """Emulate Armadillo's non-BLAS ``direct_dot_generic(x, ones)``.

    This is a diagnostic reduction, not the path used by the installed R
    build for the SCT fixtures: Armadillo delegates vectors longer than 32
    elements to BLAS there.  Keeping the implementation named and isolated
    prevents an experimental summation order from silently changing the
    stable baseline.
    """
    first = 0.0
    second = 0.0
    iterator = iter(values)
    while True:
        try:
            first += float(next(iterator))
        except StopIteration:
            break
        try:
            second += float(next(iterator))
        except StopIteration:
            break
    return first + second


def _resolve_cr_sum(cr_sum_function):
    if cr_sum_function is None or cr_sum_function == "sequential":
        return _sequential_sum
    if cr_sum_function == "armadillo_generic":
        return _armadillo_generic_dot_ones
    if callable(cr_sum_function):
        return cr_sum_function
    raise ValueError(f"unknown Cox-Reid reduction: {cr_sum_function!r}")


def _r_count_mean_variance(y: np.ndarray) -> tuple[float, float]:
    """R 4.5/MinGW-compatible mean and sample variance for UMI counts.

    The native routine preserves R's x87 ``LDOUBLE`` accumulation and the
    intermediate binary64 store used by the variance pass.
    """
    from .port_native import count_mean_variance_native

    return count_mean_variance_native(y)


def _port_scaled_newton_step(gradient: float, hessian: float,
                             scale: float = 1.0) -> tuple[float, float]:
    """One-dimensional positive-definite ``DG7QTS`` path.

    R's ``nlminb`` passes ``scale = 1`` and initializes PORT with
    ``IV(DTYPE) = 0``.  Consequently ``DRMNHB`` deliberately skips
    ``DD7DUP`` and keeps that scale fixed; it does *not* replace it with
    ``sqrt(abs(hessian))``.  The two triangular divisions below reproduce
    ``DL7IVM`` followed by ``DL7ITV``.  Collapsing them to ``-g/h`` loses a
    trajectory-relevant ULP for some genes.

    ``scale`` remains explicit so this primitive can also model a non-default
    user supplied ``nlminb(scale=...)`` value, but it is never adapted here.
    """
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("PORT scale must be finite and positive")
    reciprocal_scale = 1.0/scale
    scaled_hessian = reciprocal_scale*hessian/scale
    cholesky = math.sqrt(scaled_hessian)
    q = (gradient/scale)/cholesky
    predicted = 0.5*(q*q)
    q = q/cholesky
    return -q/scale, predicted


def _digamma_terms(y, inv, table, digamma_function=None):
    digamma_function = digamma if digamma_function is None else digamma_function
    if table is None:
        iterator = ((float(v), 1.0) for v in y)
    else:
        iterator = zip(table[0], table[1])
    dig = 0.0
    sum_y = 0.0
    sum_prod = 0.0
    max_y = 0.0
    for value, frequency in iterator:
        dig += frequency * float(digamma_function(value + inv))
        sum_y += frequency * value
        sum_prod += frequency * (value - 1.0) * value
        max_y = max(max_y, value)
    correction = sum_prod/(2.0*inv) if inv > 1e5 else 0.0
    if max_y*1e6 < inv:
        return sum_y-correction
    dig = (dig-len(y)*float(digamma_function(inv)))*inv
    return min(dig, sum_y-correction)


def loglikelihood(y, mu, log_alpha, *, cox_reid=True, table=None,
                  lgamma_function=None, cr_sum_function=None,
                  exp_function=None):
    """glmGamPoi-like adjusted log likelihood, including its raw cancellation."""
    lgamma_function = gammaln if lgamma_function is None else lgamma_function
    exp_function = math.exp if exp_function is None else exp_function
    alpha = exp_function(float(log_alpha))
    # R_pow_di(alpha, -1) first forms 1/alpha; spell that order explicitly.
    inv = 1.0/alpha
    cr = 0.0
    if cox_reid:
        weights = 1.0/(1.0/mu+alpha)
        cr_sum = _resolve_cr_sum(cr_sum_function)
        cr = -0.5*math.log(max(cr_sum(weights), 1e-50))*_CR
    if table is None:
        gamma = _sequential_sum(lgamma_function(y+inv))
    else:
        gamma = _sequential_sum(freq*float(lgamma_function(value+inv))
                                for value, freq in zip(table[0], table[1]))
    gamma -= len(y)*float(lgamma_function(inv))
    part = _sequential_sum((-float(value)-inv)*math.log(float(mean)+inv)
                           for value, mean in zip(y, mu))
    part -= len(y)*inv*math.log(alpha)
    return gamma+part+cr


def score(y, mu, log_alpha, *, cox_reid=True, table=None,
          digamma_function=None, cr_sum_function=None, exp_function=None):
    """Derivative of adjusted log likelihood with respect to log(alpha)."""
    exp_function = math.exp if exp_function is None else exp_function
    alpha = exp_function(float(log_alpha))
    inv = 1.0/alpha
    cr = 0.0
    if cox_reid:
        weights = 1.0/(1.0/mu+alpha)
        cr_sum = _resolve_cr_sum(cr_sum_function)
        b = cr_sum(weights)
        db = -cr_sum(value*value for value in weights)
        b_inv = 1.0/(b+1e-6)
        # Preserve inv_sympd(b + eps) * db before applying -1/2.  Writing the
        # algebraically reduced numerator/denominator form rounds differently.
        cr = -0.5*(b_inv*db)*_CR
    dig = _digamma_terms(y, inv, table, digamma_function)
    part = 0.0
    for value, mean in zip(y, mu):
        value, mean = float(value), float(mean)
        x = mean*alpha
        if x < 1e-10:
            term = x*x*(1.0/(1.0+x)-0.5)
        elif x < 1e-4:
            inv_x = 1.0/(1.0+x)
            upper = x*x*inv_x
            lower = x*x*(inv_x-0.5)
            suggested = math.log(1.0+x)-mean/(mean+inv)
            term = max(min(suggested, upper), lower)
        else:
            term = math.log(1.0+x)-mean/(mean+inv)
        # Upstream performs two ``+=`` statements. Combining their right-hand
        # sides changes the accumulator rounding on flat likelihoods.
        part += term
        part += value/(mean+inv)
    return part*inv-dig+cr*alpha


def score_derivative(y, mu, log_alpha, *, cox_reid=True, table=None,
                     digamma_function=None, trigamma_function=None,
                     cr_sum_function=None, exp_function=None):
    """Derivative of ``score`` with respect to log(alpha)."""
    exp_function = math.exp if exp_function is None else exp_function
    alpha = exp_function(float(log_alpha))
    # The pinned Windows R 4.5.2 binary's R_pow_di(alpha, -2) rounds like
    # 1/(alpha*alpha), despite the source-level reciprocal-then-square loop.
    # The distinction is visible after trigamma cancellation.  This is pinned
    # binary behavior, not a claim about every R compiler/toolchain.
    inv = 1.0/alpha
    inv2 = 1.0/(alpha*alpha)
    cr = 0.0
    cr2 = 0.0
    if cox_reid:
        weights = 1.0/(1.0/mu+alpha)
        cr_sum = _resolve_cr_sum(cr_sum_function)
        b = cr_sum(weights)+1e-6
        b_inv = 1.0/b
        db = -cr_sum(value*value for value in weights)
        d2b = 2.0*cr_sum(value*value*value for value in weights)
        ddet = b_inv*db
        d_i_db_square = ddet*ddet
        # glmGamPoi materializes d_i_db, then multiplies that 1x1 matrix by
        # itself.  Recomputing the quotient and applying Python ``**2`` leaves
        # a cancellation residue that is absent in the RcppArmadillo path.
        d2det = ddet*ddet-d_i_db_square+b_inv*d2b
        cr = (0.5*ddet*ddet-0.5*d2det)*_CR
        cr2 = -0.5*ddet*_CR
    digamma_function = digamma if digamma_function is None else digamma_function
    trigamma_function = (lambda value: polygamma(1, value)) if trigamma_function is None else trigamma_function
    if table is None:
        iterator = ((float(v), 1.0) for v in y)
    else:
        iterator = zip(table[0], table[1])
    dig = 0.0
    tri = 0.0
    for value, frequency in iterator:
        dig += frequency*float(digamma_function(value+inv))
        tri += frequency*float(trigamma_function(value+inv))
    dig -= len(y)*float(digamma_function(inv))
    # Keep glmGamPoi's two-statement rounding order.  The subtraction is
    # strongly cancellation-sensitive on nearly Poisson genes, so factoring
    # out inv2 changes the optimizer trajectory despite algebraic equality.
    tri *= inv2
    tri -= inv2*len(y)*float(trigamma_function(inv))
    part1 = 0.0
    part2 = 0.0
    for value, mean in zip(y, mu):
        value, mean = float(value), float(mean)
        one = 1.0+mean*alpha
        part1 += math.log(one)+(value-mean)/(mean+inv)
        # C++ source uses two left-associative divisions, not one division by
        # a squared denominator.
        part2 += (mean*mean*alpha+value)/one/one
    combined = -2.0*inv*(part1-dig)+(part2+tri)
    return combined+cr*(alpha*alpha)+(part1-dig)*inv+cr2*alpha


def fit_dispersion_newton(y, mu, *, max_iter=200, relative_tolerance=1e-10,
                          return_diagnostics=False, use_r45_mathlib=False,
                          cr_sum_function=None):
    """Isolated Newton candidate for the primary glmGamPoi ``nlminb`` path.

    It fails closed on non-Newton or unstable paths rather than silently
    claiming PORT equivalence.  The return value is alpha, not SCT theta.
    """
    y = np.asarray(y, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    if y.ndim != 1 or mu.shape != y.shape or len(y) < 2:
        raise ValueError('expected equal one-dimensional y and mu')
    if np.any(y < 0) or np.any(mu < 0) or not np.isfinite(y).all() or not np.isfinite(mu).all():
        raise ValueError('y and mu must be finite and non-negative')
    if np.all(y == 0):
        result = (0.0, {'branch': 'all_zero', 'iterations': 0, 'converged': True})
        return result if return_diagnostics else result[0]
    mu = np.where(mu == 0, 1e-6, mu)
    table = _table_if_small(y)
    if use_r45_mathlib:
        from .r45_mathlib import (digamma as digamma_function,
                                  trigamma as trigamma_function)
        if os.name == "nt":
            from .r45_mathlib import exp as exp_function
        else:
            # Linux R and glmGamPoi use the platform libm exponential.  The
            # Windows compatibility exponential is intentionally Rtools-only.
            exp_function = math.exp
        # On the pinned Windows build SciPy gammaln, not the literal Python
        # lgamma source translation, matches initialized R 4.5.2 at every
        # traced objective callback (81/81 bitwise). Keep the default here.
        lgamma_function = None
    else:
        digamma_function = trigamma_function = lgamma_function = exp_function = None
    score_options = {'table': table, 'digamma_function': digamma_function,
                     'cr_sum_function': cr_sum_function,
                     'exp_function': exp_function}
    derivative_options = {**score_options, 'trigamma_function': trigamma_function}
    if score(y, mu, math.log(1e-8), **score_options) < 0:
        result = (0.0, {'branch': 'left_score', 'iterations': 0, 'converged': True})
        return result if return_diagnostics else result[0]
    average, centered = _r_count_mean_variance(y)
    start = (centered-average)/(average*average)
    if not math.isfinite(start) or start <= 0:
        start = 0.5
    x = min(max(math.log(start), _LOWER), _UPPER)
    objective = -loglikelihood(y, mu, x, table=table, lgamma_function=lgamma_function,
                               cr_sum_function=cr_sum_function,
                               exp_function=exp_function)
    diagnostics = {'branch': 'primary_newton', 'iterations': 0, 'converged': False,
                   'message': 'iteration limit'}
    # ``nlminb`` defaults to scale=1 and PORT's DTYPE=0 keeps it fixed.
    port_scale = 1.0
    for iteration in range(1, max_iter+1):
        gradient = -score(y, mu, x, **score_options)
        hessian = -score_derivative(y, mu, x, **derivative_options)
        if not math.isfinite(gradient) or not math.isfinite(hessian) or hessian <= 0:
            diagnostics.update(iterations=iteration-1, message='non-positive or non-finite Hessian')
            break
        step, predicted = _port_scaled_newton_step(
            gradient, hessian, port_scale)
        candidate = min(max(x+step, _LOWER), _UPPER)
        candidate_objective = -loglikelihood(y, mu, candidate, table=table,
                                             lgamma_function=lgamma_function,
                                             cr_sum_function=cr_sum_function,
                                             exp_function=exp_function)
        if not math.isfinite(candidate_objective):
            diagnostics.update(iterations=iteration-1, message='non-finite objective')
            break
        old_objective = objective
        actual = old_objective-candidate_objective
        diagnostics['iterations'] = iteration
        if candidate_objective >= old_objective:
            # DA7SST restores the current point after a losing full-Newton
            # trial.  All traced flat PBMC fits then satisfy relative function
            # convergence at the restored point.
            if predicted <= relative_tolerance*abs(old_objective):
                diagnostics.update(converged=True,
                                   message='relative convergence after restored trial')
            else:
                diagnostics.update(message='losing Newton trial')
            break
        x, objective = candidate, candidate_objective
        # DA7SST assesses convergence *after* evaluating and accepting the
        # candidate. NREDUC is the full Newton model reduction at the previous
        # point, and HALF*FDIF <= PREDUC rejects misleadingly small model
        # reductions.  A pre-step test stops one iterate too early.
        if (0.5*actual <= predicted and
                predicted <= relative_tolerance*abs(old_objective)):
            diagnostics.update(converged=True, message='relative convergence candidate')
            break
    alpha = (math.exp if exp_function is None else exp_function)(x)
    result = (alpha, diagnostics)
    return result if return_diagnostics else result[0]


def fit_dispersion_trust_region(y, mu, *, max_iter=200, return_diagnostics=False,
                                use_r45_mathlib=False, cr_sum_function=None):
    """One-dimensional PORT-inspired bounded trust-region experiment.

    This implements the parts of R 4.5 ``DRMNHB``/``DA7SST`` needed to test
    the observed one-parameter path: bounded quadratic steps, radius changes,
    and the documented convergence quantities. It is not asserted to be a
    complete PORT translation and remains outside the SCT runtime.
    """
    y=np.asarray(y,dtype=np.float64);mu=np.asarray(mu,dtype=np.float64)
    if y.ndim!=1 or mu.shape!=y.shape or len(y)<2:raise ValueError('expected equal one-dimensional y and mu')
    if np.any(y<0) or np.any(mu<0) or not np.isfinite(y).all() or not np.isfinite(mu).all():raise ValueError('invalid y or mu')
    if np.all(y==0):
        result=(0.0,{'branch':'all_zero','iterations':0,'converged':True,'message':'all zero'})
        return result if return_diagnostics else result[0]
    mu=np.where(mu==0,1e-6,mu);table=_table_if_small(y)
    if use_r45_mathlib:
        from .r45_mathlib import (digamma as digamma_function,
                                  trigamma as trigamma_function)
        if os.name == "nt":
            from .r45_mathlib import exp as exp_function
        else:
            # Linux glmGamPoi resolves scalar exp through the platform libm.
            exp_function = math.exp
        lgamma_function = None
    else:
        digamma_function = trigamma_function = lgamma_function = exp_function = None
    score_options={'table':table,'digamma_function':digamma_function,
                   'cr_sum_function':cr_sum_function,
                   'exp_function':exp_function}
    derivative_options={**score_options,'trigamma_function':trigamma_function}
    if score(y,mu,math.log(1e-8),**score_options)<0:
        result=(0.0,{'branch':'left_score','iterations':0,'converged':True,'message':'left score'})
        return result if return_diagnostics else result[0]
    average,variance=_r_count_mean_variance(y)
    start=(variance-average)/(average*average);start=.5 if not math.isfinite(start) or start<=0 else start
    x=min(max(math.log(start),_LOWER),_UPPER);f=-loglikelihood(
        y,mu,x,table=table,lgamma_function=lgamma_function,
        cr_sum_function=cr_sum_function,exp_function=exp_function)
    radius=1.0/1.1;iterations=0;evaluations=1;message='iteration limit';converged=False
    # R's nlminb(scale=1) uses IV(DTYPE)=0, so this remains fixed.
    port_scale=1.0
    trajectory=[] if return_diagnostics else None
    for iteration in range(1,max_iter+1):
        iterations=iteration;g=-score(y,mu,x,**score_options);h=-score_derivative(y,mu,x,**derivative_options)
        if not math.isfinite(g) or not math.isfinite(h):message='non-finite derivative';break
        if h>0:
            full_step,full_reduction=_port_scaled_newton_step(g,h,port_scale)
        else:
            full_step=None;full_reduction=-1.0
        accepted=False
        for _attempt in range(80):
            if h>0:
                step=min(max(full_step,-radius),radius)
                is_newton=(step==full_step)
            else:
                step=-math.copysign(radius,g) if g else radius
                is_newton=False
            step=min(max(step,_LOWER-x),_UPPER-x)
            predicted=-(g*step+.5*h*step*step)
            if not math.isfinite(predicted) or predicted<=0:
                radius*=.5
                continue
            candidate=x+step;candidate_f=-loglikelihood(
                y,mu,candidate,table=table,lgamma_function=lgamma_function,
                cr_sum_function=cr_sum_function,
                exp_function=exp_function);evaluations+=1
            actual=f-candidate_f
            if trajectory is not None:
                trajectory.append({'iteration':iteration,'x':x,'objective':f,
                                   'gradient':g,'hessian':h,'scale':port_scale,
                                   'step':step,'candidate':candidate,
                                   'candidate_objective':candidate_f,
                                   'predicted':predicted,'actual':actual})
            relative_step=abs(step)/max(abs(candidate)+abs(x),1.0)
            if actual>1e-4*predicted:
                old_f=f;x,f=candidate,candidate_f;accepted=True
                ratio=actual/predicted
                if ratio<.1:radius=max(.1*radius,.5*abs(step))
                elif ratio>=.75 and not is_newton:radius=min(4.0*radius,max(2.0*abs(step),radius))
                if (full_reduction >= 0 and full_reduction <= 1e-10*abs(old_f)
                        and .5*actual <= predicted):
                    converged=True;message='relative quadratic reduction candidate'
                elif is_newton and relative_step<=math.sqrt(np.finfo(float).eps):
                    converged=True;message='x convergence candidate'
                break
            if full_reduction >= 0 and full_reduction <= 1e-10*abs(f):
                # DA7SST restores the current point when the trial full-Newton
                # step is numerically worse but its model reduction is already
                # below RFCTOL.  Do not chase floating-point objective noise by
                # repeatedly shrinking the radius.
                converged=True;message='relative convergence after rejected trial';break
            radius*=.5
            if relative_step<=100*np.finfo(float).eps:
                message='false convergence candidate';break
        if converged or not accepted:
            break
        radius=max(radius,1e-16)
    info={'branch':'primary_trust_region','iterations':iterations,'evaluations':evaluations,
          'converged':converged,'message':message,'radius':radius}
    if trajectory is not None:info['trajectory']=trajectory
    result=((math.exp if exp_function is None else exp_function)(x),info)
    return result if return_diagnostics else result[0]


def fit_dispersion_port_bfgs(y, mu, *, max_iter=200,
                             return_diagnostics=False,
                             use_r45_mathlib=False,
                             cr_sum_function=None,
                             cox_reid=True,
                             start_alpha=None):
    """Isolated gradient-only one-dimensional PORT fallback candidate.

    This mirrors the second and third ``glmGamPoi::overdispersion_mle``
    attempts.  It remains experimental until the callback-trace and fixed
    mean stress gates pass.
    """
    from .port1d import minimize_port_bfgs_1d

    y = np.asarray(y, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    if y.ndim != 1 or mu.shape != y.shape or len(y) < 2:
        raise ValueError("expected equal one-dimensional y and mu")
    if (np.any(y < 0) or np.any(mu < 0) or not np.isfinite(y).all()
            or not np.isfinite(mu).all()):
        raise ValueError("invalid y or mu")
    if np.all(y == 0):
        value = (0.0, {"branch": "all_zero", "converged": True})
        return value if return_diagnostics else value[0]
    mu = np.where(mu == 0, 1e-6, mu)
    table = _table_if_small(y)
    if use_r45_mathlib:
        from .r45_mathlib import (digamma as digamma_function,
                                  exp as exp_function)
        lgamma_function = None
    else:
        digamma_function = lgamma_function = exp_function = None
    score_options = {
        "cox_reid": cox_reid,
        "table": table,
        "digamma_function": digamma_function,
        "cr_sum_function": cr_sum_function,
        "exp_function": exp_function,
    }
    objective_options = {
        "cox_reid": cox_reid,
        "table": table,
        "lgamma_function": lgamma_function,
        "cr_sum_function": cr_sum_function,
        "exp_function": exp_function,
    }
    average, variance = _r_count_mean_variance(y)
    start = ((variance - average) / (average * average)
             if start_alpha is None else float(start_alpha))
    if not math.isfinite(start) or start <= 0.0:
        start = 0.5
    result, trajectory = minimize_port_bfgs_1d(
        lambda value: -loglikelihood(y, mu, value, **objective_options),
        lambda value: -score(y, mu, value, **score_options),
        math.log(start),
        _LOWER,
        _UPPER,
        max_iter=max_iter,
        trace=return_diagnostics,
    )
    alpha = (math.exp if exp_function is None else exp_function)(result.x)
    diagnostics = {
        "branch": "port_bfgs_cox_reid" if cox_reid else "port_bfgs_no_cox_reid",
        "iterations": result.iterations,
        "evaluations": result.evaluations,
        "converged": result.convergence == 0,
        "convergence": result.convergence,
        "message": result.message,
        "log_alpha": result.x,
    }
    if return_diagnostics:
        diagnostics["trajectory"] = trajectory
    value = (alpha, diagnostics)
    return value if return_diagnostics else value[0]


def fit_dispersion_port_native(y, mu, *, max_iter=200,
                               return_diagnostics=False,
                               use_r45_mathlib=True,
                               cr_sum_function=None):
    """Run the complete pinned glmGamPoi/R ``nlminb`` attempt chain.

    This is the strict candidate: R 4.5.2's original PORT Fortran routines
    choose and accept steps, while Python evaluates the audited glmGamPoi
    objective, score, and analytic score derivative.  A non-converged
    analytic-Hessian attempt is retried with PORT's BFGS approximation and,
    if necessary, once more without the Cox--Reid adjustment, matching
    ``conventional_overdispersion_mle``.
    """
    from .port_native import (glmgampoi_loglikelihood_native,
                              glmgampoi_score_derivative_native,
                              glmgampoi_score_native,
                              minimize_port_native_1d,
                              native_library_path)

    y = np.asarray(y, dtype=np.float64)
    mu = np.asarray(mu, dtype=np.float64)
    if y.ndim != 1 or mu.shape != y.shape or len(y) < 2:
        raise ValueError("expected equal one-dimensional y and mu")
    if (np.any(y < 0) or np.any(mu < 0) or not np.isfinite(y).all()
            or not np.isfinite(mu).all()):
        raise ValueError("invalid y or mu")
    if np.all(y == 0):
        value = (0.0, {
            "branch": "all_zero",
            "converged": True,
            "convergence": 0,
            "iterations": 0,
            "message": "All counts y are 0.",
            "stages": [],
        })
        return value if return_diagnostics else value[0]

    mu = np.where(mu == 0, 1e-6, mu)
    table = _table_if_small(y)
    if use_r45_mathlib:
        from .r45_mathlib import (digamma as digamma_function,
                                  trigamma as trigamma_function)
        if os.name == "nt":
            from .r45_mathlib import exp as exp_function
        else:
            # Linux glmGamPoi resolves scalar exp through the platform libm.
            exp_function = math.exp
        lgamma_function = None
    else:
        digamma_function = trigamma_function = lgamma_function = None
        exp_function = None

    def options(cox_reid):
        objective_options = {
            "cox_reid": cox_reid,
            "table": table,
            "lgamma_function": lgamma_function,
            "cr_sum_function": cr_sum_function,
            "exp_function": exp_function,
        }
        score_options = {
            "cox_reid": cox_reid,
            "table": table,
            "digamma_function": digamma_function,
            "cr_sum_function": cr_sum_function,
            "exp_function": exp_function,
        }
        derivative_options = {
            **score_options,
            "trigamma_function": trigamma_function,
        }
        return objective_options, score_options, derivative_options

    _, left_score_options, _ = options(True)
    if score(y, mu, math.log(1e-8), **left_score_options) < 0:
        value = (0.0, {
            "branch": "left_score",
            "converged": True,
            "convergence": 0,
            "iterations": 0,
            "message": "Even for very small theta, no maximum identified",
            "stages": [],
        })
        return value if return_diagnostics else value[0]

    from .port_native import count_dispersion_start_native
    start_alpha, start = count_dispersion_start_native(y)
    start = min(max(start, _LOWER), _UPPER)

    stages = []
    callback_trace = [] if return_diagnostics else None
    final = None
    attempts = (
        ("analytic_hessian_cox_reid", True, True),
        ("bfgs_cox_reid", True, False),
        ("bfgs_no_cox_reid", False, False),
    )
    for name, cox_reid, analytic_hessian in attempts:
        objective_options, score_options, derivative_options = options(cox_reid)
        def objective(value, opts=objective_options, stage=name, cr=cox_reid):
            if os.name == "nt":
                result = -loglikelihood(y, mu, value, **opts)
            else:
                result = -glmgampoi_loglikelihood_native(
                    y, mu, value, cox_reid=cr, table=table)
            if callback_trace is not None:
                callback_trace.append((stage, "objective", value, result))
            return result

        def gradient(value, cr=cox_reid, stage=name):
            result = -glmgampoi_score_native(
                y, mu,
                (math.exp if exp_function is None else exp_function)(value),
                cox_reid=cr, table=table)
            if callback_trace is not None:
                callback_trace.append((stage, "gradient", value, result))
            return result
        hessian = None
        if analytic_hessian:
            def hessian(value, cr=cox_reid, opts=derivative_options,
                        stage=name):
                result = -glmgampoi_score_derivative_native(
                    y,
                    mu,
                    (math.exp if exp_function is None else exp_function)(
                        value),
                    cox_reid=cr,
                    table=table,
                )
                if callback_trace is not None:
                    callback_trace.append((stage, "hessian", value, result))
                return result
        result = minimize_port_native_1d(
            objective,
            gradient,
            start,
            _LOWER,
            _UPPER,
            hessian=hessian,
            max_iter=max_iter,
        )
        stages.append({
            "name": name,
            "status": result.status,
            "convergence": result.convergence,
            "message": result.message,
            "iterations": result.iterations,
            "function_evaluations": result.function_evaluations,
            "gradient_evaluations": result.gradient_evaluations,
            "log_alpha": result.x,
        })
        final = result
        if result.convergence == 0:
            break

    assert final is not None
    alpha = (math.exp if exp_function is None else exp_function)(final.x)
    diagnostics = {
        "branch": stages[-1]["name"],
        "converged": final.convergence == 0,
        "convergence": final.convergence,
        "iterations": final.iterations,
        "message": final.message,
        "start_alpha": start_alpha,
        "log_alpha": final.x,
        "native_library": str(native_library_path()),
        "stages": stages,
    }
    if callback_trace is not None:
        diagnostics["callback_trace"] = callback_trace
    value = (alpha, diagnostics)
    return value if return_diagnostics else value[0]
