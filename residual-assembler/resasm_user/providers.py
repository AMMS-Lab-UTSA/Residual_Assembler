"""Residual + tangent providers for the three user paths.

Path A  python      -> a residual(u, params, state, time) callable (+ optional
                       tangent(...)); the framework runs OTILib internally.
Path C  executable  -> a black-box command that returns R^(p) in a response file;
                       the framework only solves T U^(p) = -R^(p).

Tangent sources: python function, a saved file (.npz/.npy), or the black-box
response. Everything runs locally on the user's machine.
"""

from __future__ import annotations

import importlib.util
import json
import os
import shlex
import subprocess
import tempfile
from typing import Any, Callable, Dict, List, Optional

import numpy as np

from .config import UserConfig, ConfigError


# --------------------------------------------------------------------------- #
def _load_callable(module_path: str, func_name: str) -> Callable:
    if not os.path.exists(module_path):
        raise ConfigError("residual module not found: %s" % module_path,
                          "residual:\n  module: user_residual.py")
    spec = importlib.util.spec_from_file_location("user_residual_mod", module_path)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)                # type: ignore[union-attr]
    except Exception as exc:                        # noqa: BLE001
        raise ConfigError("could not import %s: %s" % (module_path, exc))
    fn = getattr(mod, func_name, None)
    if not callable(fn):
        raise ConfigError(
            "function '%s' not found in %s" % (func_name, os.path.basename(module_path)),
            "def %s(u, params, state=None, time=None):\n    ...\n    return R" % func_name)
    return fn


class PythonResidual:
    """Wraps a user residual callable (and optional tangent callable)."""

    kind = "python"

    def __init__(self, cfg: UserConfig):
        mod = cfg.path(cfg.residual["module"])
        self.residual = _load_callable(mod, cfg.residual.get("function", "residual"))
        self.tangent_fn: Optional[Callable] = None
        tan = cfg.tangent or {}
        if str(tan.get("type", "")).lower() == "python":
            fname = tan.get("function", "tangent")
            self.tangent_fn = _load_callable(mod, fname)

    def eval_residual(self, u, params, state=None, time=None):
        return self.residual(u, params, state, time)

    def eval_tangent(self, u, params, state=None, time=None):
        if self.tangent_fn is None:
            return None
        return np.asarray(self.tangent_fn(u, params, state, time), float)


class ExecutableResidual:
    """Black-box: a command that reads request.json and writes response.npz.

    Request JSON fields: schema, u, parameters, seed_directions, basis_count,
    truncation_order, order, direction_map, u_star_coefficients, state, time,
    dtime. Response (.npz or .json): R_order_<p> arrays, optional `tangent`,
    optional `diagnostics` (JSON string)."""

    kind = "executable"

    def __init__(self, cfg: UserConfig):
        self.cfg = cfg
        self.command = cfg.residual["command"]

    def eval_rhs(self, u, parameters: Dict[str, float], order: int,
                 direction_map: Dict[int, List[List[int]]],
                 u_star_coefficients: Optional[Dict[int, Any]] = None,
                 state=None, time=(0.0, 0.0), dtime=0.0) -> Dict[str, Any]:
        names = list(parameters.keys())
        request = {
            "schema": "resasm-user-request/1",
            "u": [float(x) for x in np.asarray(u, float).ravel()],
            "parameters": {k: float(v) for k, v in parameters.items()},
            "seed_directions": {n: i + 1 for i, n in enumerate(names)},
            "basis_count": len(names),
            "truncation_order": int(order),
            "order": int(order),
            "direction_map": {str(p): v for p, v in direction_map.items()},
            "u_star_coefficients": u_star_coefficients or {},
            "state": state, "time": list(time), "dtime": float(dtime),
        }
        workdir = self.cfg.base_dir
        with tempfile.TemporaryDirectory(dir=workdir) as td:
            req = os.path.join(td, "request.json")
            resp = os.path.join(td, "response.npz")
            with open(req, "w", encoding="utf-8") as fh:
                json.dump(request, fh)
            argv = self._build_argv(req, resp)
            proc = subprocess.run(argv, cwd=workdir, capture_output=True,
                                  text=True, timeout=3600)
            if not (os.path.exists(resp) or os.path.exists(resp[:-4] + ".json")):
                raise ConfigError(
                    "black-box executable produced no response file.\n"
                    "command: %s\nstdout:\n%s\nstderr:\n%s"
                    % (self.command, proc.stdout[-1500:], proc.stderr[-1500:]))
            return _read_response(resp)

    def _build_argv(self, req: str, resp: str) -> List[str]:
        cmd = self.command.replace("{request}", req).replace("{response}", resp)
        parts = shlex.split(cmd, posix=(os.name != "nt"))
        # route a .bat launcher through cmd /c on Windows
        if os.name == "nt" and parts and parts[0].lower().endswith(".bat"):
            return ["cmd", "/c"] + parts
        return parts


def _read_response(resp_npz: str) -> Dict[str, Any]:
    out: Dict[str, Any] = {"orders": {}, "tangent": None, "diagnostics": {}}
    if os.path.exists(resp_npz):
        data = np.load(resp_npz, allow_pickle=True)
        for key in data.files:
            if key.startswith("R_order_"):
                p = int(key.rsplit("_", 1)[1])
                out["orders"][p] = np.asarray(data[key], float)
            elif key == "tangent":
                out["tangent"] = np.asarray(data[key], float)
            elif key == "diagnostics":
                try:
                    out["diagnostics"] = json.loads(str(data[key]))
                except Exception:                   # noqa: BLE001
                    out["diagnostics"] = {"raw": str(data[key])}
        return out
    resp_json = resp_npz[:-4] + ".json"
    with open(resp_json, "r", encoding="utf-8") as fh:
        d = json.load(fh)
    rbo = d.get("residual_coefficients_by_order", {})
    for p, arr in rbo.items():
        out["orders"][int(p)] = np.asarray(arr, float)
    if d.get("tangent") is not None:
        out["tangent"] = np.asarray(d["tangent"], float)
    out["diagnostics"] = d.get("diagnostics", {})
    return out


def load_tangent_file(path: str) -> np.ndarray:
    if path.endswith(".npz"):
        data = np.load(path, allow_pickle=True)
        key = "tangent" if "tangent" in data.files else data.files[0]
        return np.asarray(data[key], float)
    return np.asarray(np.load(path, allow_pickle=True), float)
