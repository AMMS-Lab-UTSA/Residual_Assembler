"""Orchestration: run a complete sensitivity job from a resasm.yml."""

from __future__ import annotations

import os
import time as _time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np

from .config import UserConfig, ConfigError, load_config
from .providers import (PythonResidual, ExecutableResidual, load_tangent_file)
from . import oti_global


@dataclass
class RunResult:
    ok: bool
    output_dir: str
    private_dir: str
    public_dir: str
    summary: Dict[str, Any]


# --------------------------------------------------------------------------- #
def _load_solution(cfg: UserConfig) -> np.ndarray:
    sol = cfg.solution
    if sol.get("values") is not None:
        u = np.asarray(sol["values"], float).ravel()
    else:
        p = cfg.path(sol["file"])
        if not os.path.exists(p):
            raise ConfigError("solution file not found: %s" % p,
                              "solution:\n  file: solution.npy")
        u = np.asarray(np.load(p, allow_pickle=True), float).ravel()
    if u.size != cfg.unknowns:
        raise ConfigError(
            "solution vector has %d entries but problem.unknowns is %d"
            % (u.size, cfg.unknowns),
            "problem:\n  unknowns: %d" % u.size)
    return u


def _load_state(cfg: UserConfig):
    st = cfg.state or {}
    if st.get("file"):
        return np.load(cfg.path(st["file"]), allow_pickle=True)
    return st.get("values")


def _time_pair(cfg: UserConfig):
    t = cfg.time or {}
    return (float(t.get("time", 0.0)), float(t.get("time", 0.0))), float(t.get("dtime", 0.0))


def _free_mask(cfg: UserConfig, ndof: int) -> np.ndarray:
    c = cfg.constraints or {}
    free = np.ones(ndof, bool)
    if c.get("prescribed"):
        for i in c["prescribed"]:
            free[int(i)] = False
    elif c.get("free"):
        free[:] = False
        for i in c["free"]:
            free[int(i)] = True
    return free


def _acquire_tangent(cfg, provider, u, params, state, time_pair, ndof):
    tan = cfg.tangent or {}
    ttype = str(tan.get("type", "")).lower()
    if ttype == "file" and tan.get("file"):
        T = load_tangent_file(cfg.path(tan["file"]))
        return _check_tangent(T, ndof), "file"
    if ttype == "python" and isinstance(provider, PythonResidual):
        T = provider.eval_tangent(u, params, state, time_pair[0])
        if T is not None:
            return _check_tangent(T, ndof), "python"
    return None, ttype or "response"


def _check_tangent(T, ndof):
    T = np.asarray(T, float)
    if T.shape != (ndof, ndof):
        raise ConfigError(
            "tangent has shape %s but must be (%d, %d)" % (T.shape, ndof, ndof),
            "tangent:\n  type: file\n  file: tangent.npz")
    return T


def _wants_rhs_fd(validation: Optional[Dict[str, Any]]) -> bool:
    """True if the user asked for the RHS finite-difference cross-check.

    Canonical key: `validation.rhs_finite_difference_check`.
    Legacy key `validation.finite_difference` is still accepted (it meant the
    same thing) but the output always reports the accurate name.
    """
    v = validation or {}
    return bool(v.get("rhs_finite_difference_check", v.get("finite_difference")))


def _solution_fd_note(validation: Optional[Dict[str, Any]]) -> Optional[str]:
    """`validation.solution_finite_difference_solver` is a RESERVED field for a
    future true solution-level FD check (re-solving the nonlinear problem at
    perturbed parameters). It is not implemented; say so plainly rather than
    silently ignoring it."""
    v = validation or {}
    if v.get("solution_finite_difference_solver"):
        return ("validation.solution_finite_difference_solver is reserved but NOT "
                "implemented: no solution-level finite-difference validation was "
                "performed. Only the RHS (residual-derivative) check ran.")
    return None


# --------------------------------------------------------------------------- #
def run_from_config(path: str) -> RunResult:
    """Run a job. Dispatches on the config kind:

      * **Path A (assembly recipe)** -- the config names a `mesh:`; we assemble R
        from the ingredients. THIS IS THE PRIMARY PATH.
      * Path B/C -- the config names a `residual:` (executable / python).
    """
    # field-driven residual sensitivity (M2) also names a `mesh:`, so it must be
    # detected BEFORE the generic assembly-recipe check below.
    from .field_recipe import is_field_recipe, run_field_recipe
    if is_field_recipe(path):
        return run_field_recipe(path)

    from .checks import _is_recipe
    if _is_recipe(path):
        from .assembly_runner import run_recipe
        return run_recipe(path)
    return _run_residual_config(path)


def _run_residual_config(path: str) -> RunResult:
    t0 = _time.time()
    cfg = load_config(path)
    u = _load_solution(cfg)
    ndof = u.size
    params = dict(cfg.parameters)
    names = list(params.keys())
    state = _load_state(cfg)
    time_pair, dtime = _time_pair(cfg)
    free = _free_mask(cfg, ndof)

    # None == "the real residual was never available" (black-box). Never fabricate
    # a zero vector here: a reader of private/residual_real.npz would see ||R||=0
    # and conclude equilibrium was verified when nothing of the sort happened.
    R_real = None
    validation: Dict[str, Any] = {"status": "ok"}

    if cfg.residual_type in ("python", "element"):
        provider = PythonResidual(cfg)
        # real residual + free-DOF norm
        r = provider.eval_residual(u, params, state, time_pair[0])
        R_real = np.asarray([float(x) for x in r], float)
        validation["residual_free_norm"] = float(np.linalg.norm(R_real[free]))
        T, tsrc = _acquire_tangent(cfg, provider, u, params, state, time_pair, ndof)
        if T is None:
            raise ConfigError(
                "no tangent available for the Python path",
                "tangent:\n  type: python\n  function: tangent\n"
                "# or:\ntangent:\n  type: file\n  file: tangent.npz")
        R_orders, U_orders, diag = oti_global.solve_python(
            provider.eval_residual, u, params, cfg.order, T, free,
            state, time_pair, dtime)
        tangent_used = T
        if _wants_rhs_fd(cfg.validation):
            validation["rhs_finite_difference_check"] = rhs_finite_difference_check(
                provider, u, params, T, free, state, time_pair, U_orders.get(1))
    elif cfg.residual_type == "executable":
        exe = ExecutableResidual(cfg)
        T, tsrc = _acquire_tangent(cfg, None, u, params, state, time_pair, ndof)
        R_orders, U_orders, tangent_used, diag = oti_global.solve_executable(
            exe, u, params, cfg.order, T, free, state, time_pair, dtime)
        tsrc = tsrc if T is not None else "response"
        # A black-box solver returns only the perturbed RHS; it never hands back
        # the real residual, so there is NOTHING to take a norm of. Keep this
        # None and say why -- never report a fabricated 0.0 (that would imply we
        # verified equilibrium when we did not).
        validation["residual_free_norm"] = None
        validation["residual_free_norm_reason"] = (
            "black-box did not expose real residual")
    else:
        raise ConfigError("unsupported residual.type: %s" % cfg.residual_type)

    note = _solution_fd_note(cfg.validation)
    if note:
        validation.setdefault("notes", []).append(note)

    labels = oti_global.direction_labels(len(names), cfg.order)
    result = {
        "parameter_names": names, "ndof": ndof, "free_mask": free,
        "R_real": R_real, "tangent": tangent_used, "tangent_source": tsrc,
        "R_orders": R_orders, "U_orders": U_orders, "direction_labels": labels,
        "validation": validation, "diagnostics": diag,
        "timing": {"total_seconds": round(_time.time() - t0, 4),
                   "orders": cfg.order, "unknowns": ndof},
    }
    from .output import write_outputs
    out_dir = cfg.path(cfg.output_dir)
    paths = write_outputs(out_dir, cfg, result)
    return RunResult(ok=True, output_dir=out_dir, private_dir=paths["private"],
                     public_dir=paths["public"],
                     summary={"parameters": names, "order": cfg.order,
                              "tangent_source": tsrc,
                              "residual_free_norm": validation.get("residual_free_norm"),
                              "orders_solved": sorted(U_orders.keys())})


def rhs_finite_difference_check(provider, u, params, T, free, state, time_pair, U1,
                                h_rel=1e-6) -> Dict[str, Any]:
    """Finite-difference check of the RHS (residual derivative) — NOT a re-solve.

    What this DOES check
    --------------------
    For each parameter a_i it central-differences the residual with respect to
    a_i **at the fixed solution u**::

        dR/da_i ~ (R(u, a_i+h) - R(u, a_i-h)) / 2h

    then solves the SAME linear system with the SAME tangent T that the
    hypercomplex path used::

        u1_fd = -T_ff^-1 (dR/da_i)_f

    and compares that to the order-1 sensitivity U^(1). So it validates the
    generated RHS / residual derivative and the solve that consumes it.

    What this does NOT check
    ------------------------
    It does NOT re-solve the nonlinear problem at the perturbed parameters, so it
    cannot detect an error in the tangent itself (a wrong T cancels out: it is
    used on both sides), nor anything about the converged solution u. It is
    therefore NOT a "full solution finite-difference validation" and must not be
    reported as one. A true solution-level FD check would require re-running the
    user's nonlinear solver at a_i +/- h (see the reserved config field
    `validation.solution_finite_difference_solver`, not implemented).
    """
    if U1 is None:
        return {"status": "skipped (no order-1 solution)",
                "checks": "d(residual)/d(parameter) at fixed u"}
    names = list(params.keys())
    Tff = np.asarray(T, float)[np.ix_(free, free)]
    ndof = u.size
    max_rel = 0.0
    for i, n in enumerate(names):
        v0 = float(params[n])
        h = h_rel * max(1.0, abs(v0))
        pp = dict(params); pp[n] = v0 + h
        pm = dict(params); pm[n] = v0 - h
        rp = np.asarray([float(x) for x in provider.eval_residual(u, pp, state, time_pair[0])], float)
        rm = np.asarray([float(x) for x in provider.eval_residual(u, pm, state, time_pair[0])], float)
        dRda = (rp - rm) / (2.0 * h)
        u1_fd = np.zeros(ndof)
        u1_fd[free] = np.linalg.solve(Tff, -dRda[free])
        denom = np.linalg.norm(U1[:, i]) or 1.0
        rel = float(np.linalg.norm(u1_fd - U1[:, i]) / denom)
        max_rel = max(max_rel, rel)
    return {"max_rel_error": max_rel,
            "status": "pass" if max_rel < 1e-4 else "warn",
            "checks": "d(residual)/d(parameter) at fixed u, solved with the same "
                      "tangent (does NOT re-solve the nonlinear problem)"}
