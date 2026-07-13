"""`resasm check` — a step-by-step readiness report that stops at the first
actionable missing item. Each line is `[ok]`, `[warn]`, or `[fail]`.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from .config import UserConfig, ConfigError, load_config
from .providers import PythonResidual, ExecutableResidual, load_tangent_file
from . import oti_global
from . import runner as _runner


@dataclass
class CheckReport:
    ok: bool = True
    lines: List[str] = field(default_factory=list)

    def render(self) -> str:
        return "\n".join(self.lines)


def _ok(rep, msg):
    rep.lines.append("[ok] " + msg)


def _warn(rep, msg):
    rep.lines.append("[warn] " + msg)


def _fail(rep, msg):
    rep.ok = False
    rep.lines.append("[fail] " + msg)


def check_config(path: str, warn_tol: float = 1e-4) -> CheckReport:
    """Validate everything needed for a run, stopping at the first blocker."""
    rep = CheckReport()
    try:
        cfg = load_config(path)
    except ConfigError as exc:
        _fail(rep, str(exc))
        return rep
    _ok(rep, "config parsed: %s (order %d, backend %s)"
        % (cfg.name, cfg.order, cfg.backend))

    # solution
    try:
        u = _runner._load_solution(cfg)
    except ConfigError as exc:
        _fail(rep, str(exc))
        return rep
    _ok(rep, "loaded solution vector: shape (%d,)" % u.size)

    # parameters
    _ok(rep, "loaded parameter map: %d parameters (%s)"
        % (len(cfg.parameters), ", ".join(cfg.parameters)))

    ndof = u.size
    state = _runner._load_state(cfg)
    time_pair, dtime = _runner._time_pair(cfg)
    free = _runner._free_mask(cfg, ndof)

    provider = None
    R_real = None
    if cfg.residual_type in ("python", "element"):
        try:
            provider = PythonResidual(cfg)
            r = provider.eval_residual(u, cfg.parameters, state, time_pair[0])
            R_real = np.asarray([float(x) for x in r], float)
        except ConfigError as exc:
            _fail(rep, str(exc))
            return rep
        except Exception as exc:                    # noqa: BLE001
            _fail(rep, "residual evaluation failed: %s\n"
                       "Check that residual(u, params, state, time) returns a "
                       "length-%d vector." % (exc, ndof))
            return rep
        if R_real.size != ndof:
            _fail(rep, "residual returned shape (%d,) but problem.unknowns is %d"
                  % (R_real.size, ndof))
            return rep
        _ok(rep, "residual evaluated: shape (%d,)" % R_real.size)
    else:
        try:
            ExecutableResidual(cfg)
        except ConfigError as exc:
            _fail(rep, str(exc))
            return rep
        _ok(rep, "black-box executable configured: %s" % cfg.residual["command"])

    # tangent
    T = None
    try:
        T, tsrc = _runner._acquire_tangent(cfg, provider, u, cfg.parameters,
                                           state, time_pair, ndof)
    except ConfigError as exc:
        _fail(rep, str(exc))
        return rep
    if T is not None:
        _ok(rep, "tangent loaded (%s): shape (%d, %d)" % (tsrc, ndof, ndof))
    elif cfg.residual_type == "executable":
        _warn(rep, "tangent will come from the black-box response")
    else:
        _fail(rep, "no tangent available.\nAdd one of:\n"
                   "    tangent:\n      type: python\n      function: tangent\n"
                   "or\n    tangent:\n      type: file\n      file: tangent.npz")
        return rep

    # OTILib (only the internal python/element path needs it)
    if cfg.residual_type in ("python", "element"):
        from residual_core.algebra.otilib_adapter import (otilib_available,
                                                          otilib_status)
        if cfg.backend in ("otilib", "oti", "hypad"):
            if otilib_available():
                _ok(rep, "OTILib available: order %d, basis %d"
                    % (cfg.order, len(cfg.parameters)))
            else:
                _fail(rep, otilib_status()["error"])
                return rep
        elif cfg.backend in ("dual1", "dual") and cfg.order > 1:
            _fail(rep, "backend 'dual1' supports order 1 only; set backend: otilib "
                       "for order %d." % cfg.order)
            return rep

    # residual norm quality
    if R_real is not None:
        rn = float(np.linalg.norm(R_real[free]))
        if rn > warn_tol:
            _warn(rep, "residual norm on free DOFs is %.3e (is u converged?)" % rn)
        else:
            _ok(rep, "residual norm on free DOFs is %.3e" % rn)

    # tiny order-1 RHS + solve
    try:
        if cfg.residual_type in ("python", "element"):
            R1, U1, _d = oti_global.solve_python(
                provider.eval_residual, u, cfg.parameters, 1, T, free,
                state, time_pair, dtime)
        else:
            exe = ExecutableResidual(cfg)
            R1, U1, _T, _d = oti_global.solve_executable(
                exe, u, cfg.parameters, 1, T, free, state, time_pair, dtime)
    except ConfigError as exc:
        _fail(rep, str(exc))
        return rep
    except Exception as exc:                        # noqa: BLE001
        _fail(rep, "RHS generation / solve failed: %s" % exc)
        return rep
    m = len(cfg.parameters)
    _ok(rep, "RHS order 1 generated: shape (%d, %d)" % (ndof, m))
    _ok(rep, "sensitivity solve completed")
    return rep
