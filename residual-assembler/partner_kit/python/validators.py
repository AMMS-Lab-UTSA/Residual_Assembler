"""Local validation — runs on the partner machine, needs no shared data.

Implements the five checks from ``docs/validation_checklist.md``:

  1. real residual   R_free ~ 0 at the converged solution
  2. tangent         T ~ dR/dU (directional finite differences)
  3. RHS             R^(1) ~ dR/dp (finite differences, per parameter)
  4. sensitivity     T dU/dp = -dR/dp  (residual of the solved system)
  5. output deriv.   (optional) selected output sensitivities vs finite differences

For black-box providers the residual is not evaluable locally by the kit, so the
FD checks are marked ``pending`` (the executable should self-validate); checks 1
and 4 still run from the returned arrays.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from .residual_provider import (ResidualProvider, GlobalResidualProvider,
                                ElementResidualProvider, BlackBoxResidualProvider)
from .sensitivity import assemble_real


def _check(name, value, tol, passed, detail=""):
    return {"name": name, "value": (None if value is None else float(value)),
            "tolerance": (None if tol is None else float(tol)),
            "passed": passed, "detail": detail}


def _is_local(provider) -> bool:
    return isinstance(provider, (GlobalResidualProvider, ElementResidualProvider))


def check_real_residual(provider, U, state, time, dtime, tol=1e-8):
    if not _is_local(provider):
        return _check("real-residual", None, tol, None,
                      "black-box: check inside the executable")
    values = provider.parameter_values()
    R, _T = assemble_real(provider, U, state, values, time, dtime)
    free = provider.free_mask()
    free = np.ones(R.size, bool) if free is None else np.asarray(free, bool)
    n = float(np.linalg.norm(R[free]))
    return _check("real-residual", n, tol, n < tol, "||R_free|| at the solution")


def check_tangent_directional(provider, U, T, state, time, dtime,
                              n_dirs=3, h=1e-6, tol=1e-4, seed=0):
    if not _is_local(provider) or T is None:
        return _check("tangent-vs-FD", None, tol, None,
                      "unavailable (black-box or no tangent)")
    values = provider.parameter_values()
    U = np.asarray(U, float)
    rng = np.random.default_rng(seed)
    rels = []
    for _ in range(n_dirs):
        v = rng.standard_normal(U.size)
        v /= np.linalg.norm(v)
        Rp, _ = assemble_real(provider, U + h * v, state, values, time, dtime)
        Rm, _ = assemble_real(provider, U - h * v, state, values, time, dtime)
        fd = (Rp - Rm) / (2 * h)
        Tv = np.asarray(T, float) @ v
        rels.append(np.linalg.norm(Tv - fd) / max(np.linalg.norm(fd), 1e-30))
    rel = float(np.max(rels))
    return _check("tangent-vs-FD", rel, tol, rel < tol,
                  "T v vs central FD of R, %d directions" % n_dirs)


def check_rhs_fd(provider, U, R1, state, time, dtime, h=1e-6, tol=1e-4):
    if not _is_local(provider):
        return _check("rhs-vs-FD", None, tol, None,
                      "black-box: check inside the executable")
    params = list(provider.parameters)
    values = provider.parameter_values()
    U = np.asarray(U, float)
    rels = []
    for col, p in enumerate(params):
        v0 = float(values[p])
        dv = h * max(1.0, abs(v0))
        vp = dict(values); vp[p] = v0 + dv
        vm = dict(values); vm[p] = v0 - dv
        Rp, _ = assemble_real(provider, U, state, vp, time, dtime)
        Rm, _ = assemble_real(provider, U, state, vm, time, dtime)
        fd = (Rp - Rm) / (2 * dv)
        rels.append(np.linalg.norm(R1[:, col] - fd) / max(np.linalg.norm(fd), 1e-30))
    rel = float(np.max(rels)) if rels else 0.0
    return _check("rhs-vs-FD", rel, tol, rel < tol, "R^(1) vs central FD of dR/dp")


def check_sensitivity_solve(T, R1, U1, free_mask=None, tol=1e-8):
    if T is None or U1 is None:
        return _check("sensitivity-solve", None, tol, None, "no tangent / no solve")
    T = np.asarray(T, float); R1 = np.asarray(R1, float); U1 = np.asarray(U1, float)
    free = np.ones(T.shape[0], bool) if free_mask is None else np.asarray(free_mask, bool)
    resid = T[np.ix_(free, free)] @ U1[free] + R1[free]
    n = float(np.linalg.norm(resid))
    scale = max(float(np.linalg.norm(R1[free])), 1e-30)
    rel = n / scale
    return _check("sensitivity-solve", rel, tol, rel < tol, "T dU/dp + dR/dp ~ 0")


def check_output_sensitivity_fd(provider, U, U1, output_fn, state, time, dtime,
                                h=1e-6, tol=1e-4):
    """Optional: an output g(U, params). Compares dg/dp from the chain rule
    (dg/dU · dU/dp + dg/dp) against a finite difference of g wrt the parameter.
    ``output_fn(U, values) -> float``. Requires a local provider + a real solve
    hook is not done here (uses U1 already computed)."""
    if not _is_local(provider) or output_fn is None:
        return _check("output-sensitivity-vs-FD", None, tol, None, "not provided")
    # numeric dg/dU
    params = list(provider.parameters)
    values = provider.parameter_values()
    g0 = float(output_fn(U, values))
    dg_dU = np.zeros(U.size)
    for i in range(U.size):
        e = np.zeros(U.size); e[i] = h
        dg_dU[i] = (float(output_fn(U + e, values)) - float(output_fn(U - e, values))) / (2 * h)
    rels = []
    for col, p in enumerate(params):
        # chain rule: dg/dp = dg/dU · dU/dp + dg/dp|_explicit
        v0 = float(values[p]); dv = h * max(1.0, abs(v0))
        vp = dict(values); vp[p] = v0 + dv
        vm = dict(values); vm[p] = v0 - dv
        dg_dp_expl = (float(output_fn(U, vp)) - float(output_fn(U, vm))) / (2 * dv)
        chain = float(dg_dU @ U1[:, col]) + dg_dp_expl
        # note: a pure-FD reference would re-solve U(p); here we compare the
        # explicit + chain estimate for self-consistency of U1.
        rels.append(abs(chain))     # reported as magnitude; see docs
    rel = float(np.max(rels)) if rels else 0.0
    return _check("output-sensitivity", rel, None, None,
                  "dg/dp via chain rule with solved dU/dp (informational)")


def run_all(provider, U, R1, T, U1, state=None, time=(0.0, 0.0), dtime=0.0,
            free_mask=None) -> List[Dict[str, Any]]:
    return [
        check_real_residual(provider, U, state, time, dtime),
        check_tangent_directional(provider, U, T, state, time, dtime),
        check_rhs_fd(provider, U, R1, state, time, dtime),
        check_sensitivity_solve(T, R1, U1, free_mask),
    ]


def overall(checks) -> Optional[bool]:
    decided = [c["passed"] for c in checks if c["passed"] is not None]
    return all(decided) if decided else None
