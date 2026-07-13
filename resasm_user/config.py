"""User config (`resasm.yml`) loading + validation with actionable errors.

Minimal required shape (see docs/minimal_user_config.md):

    problem:      {name}                 # `unknowns` is OPTIONAL (inferred)
    residual:     {type: python|executable|element, ...}
    parameters:   {name: value, ...}
    solution:     {file: ...}   (or {values: [...]})
    sensitivity:  {order, backend}

`problem.unknowns` is inferred from the solution vector when omitted; if both are
given they must agree. Everything else is optional. Errors are short and tell the
user exactly what to add.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from . import _miniyaml


class ConfigError(Exception):
    """A user-facing configuration error: short message + a concrete fix."""

    def __init__(self, message: str, fix: Optional[str] = None):
        self.message = message
        self.fix = fix
        full = message
        if fix:
            full += "\n\nAdd this to resasm.yml:\n\n" + _indent(fix)
        super().__init__(full)


def _indent(block: str, n: int = 4) -> str:
    pad = " " * n
    return "\n".join(pad + ln for ln in block.strip("\n").splitlines())


@dataclass
class UserConfig:
    raw: Dict[str, Any]
    base_dir: str
    name: str
    unknowns: int
    residual_type: str
    parameters: Dict[str, float]
    order: int
    backend: str
    residual: Dict[str, Any] = field(default_factory=dict)
    tangent: Dict[str, Any] = field(default_factory=dict)
    solution: Dict[str, Any] = field(default_factory=dict)
    constraints: Dict[str, Any] = field(default_factory=dict)
    state: Dict[str, Any] = field(default_factory=dict)
    time: Dict[str, Any] = field(default_factory=dict)
    output_dir: str = "resasm_output"
    validation: Dict[str, Any] = field(default_factory=dict)

    def path(self, rel: Optional[str]) -> Optional[str]:
        if not rel:
            return None
        return rel if os.path.isabs(rel) else os.path.join(self.base_dir, rel)


def load_config(path: str) -> UserConfig:
    """Parse + validate a resasm.yml. Raises ConfigError with an actionable
    message on the first problem found."""
    if not os.path.exists(path):
        raise ConfigError(
            "Config file not found: %s" % path,
            "problem:\n  name: my_job\n  unknowns: 1\nresidual:\n  type: python\n"
            "  module: user_residual.py\n  function: residual\nparameters:\n  k: 2.0\n"
            "solution:\n  file: solution.npy\nsensitivity:\n  order: 1\n  backend: otilib")
    try:
        raw = _miniyaml.load_file(path)
    except Exception as exc:                       # noqa: BLE001
        raise ConfigError("Could not parse config %s: %s" % (path, exc))
    if not isinstance(raw, dict):
        raise ConfigError("Config must be a mapping (key: value pairs), got %s"
                          % type(raw).__name__)

    base_dir = os.path.dirname(os.path.abspath(path))

    problem = _require(raw, "problem", "problem:\n  name: my_job")
    name = _require(problem, "name", "problem:\n  name: my_job", where="problem")

    # `problem.unknowns` is OPTIONAL: it is inferred from the solution vector when
    # omitted. If BOTH are given they must agree (a silent mismatch here would
    # corrupt every downstream shape check).
    declared = problem.get("unknowns")
    if declared is not None:
        try:
            declared = int(declared)
        except (TypeError, ValueError):
            raise ConfigError("problem.unknowns must be an integer (number of DOFs)",
                              "problem:\n  unknowns: 1")
        if declared < 1:
            raise ConfigError("problem.unknowns must be >= 1 (got %d)" % declared,
                              "problem:\n  unknowns: 1")

    residual = _require(raw, "residual",
                        "residual:\n  type: python\n  module: user_residual.py\n"
                        "  function: residual")
    rtype = _require(residual, "type",
                     "residual:\n  type: python   # python | executable | element",
                     where="residual")
    rtype = str(rtype).lower()
    if rtype not in ("python", "executable", "element"):
        raise ConfigError(
            "residual.type '%s' is not supported" % rtype,
            "residual:\n  type: python   # one of: python, executable, element")
    _validate_residual(residual, rtype)

    params = _require(raw, "parameters",
                      "parameters:\n  k: 2.0\n  f: 16.0")
    if not isinstance(params, dict) or not params:
        raise ConfigError(
            "parameters must be a non-empty mapping of name -> real value",
            "parameters:\n  k: 2.0\n  f: 16.0")
    parameters = {}
    for k, v in params.items():
        try:
            parameters[str(k)] = float(v)
        except (TypeError, ValueError):
            raise ConfigError(
                "parameter '%s' must be a real number, got %r" % (k, v),
                "parameters:\n  %s: 1.0" % k)

    solution = _require(raw, "solution", "solution:\n  file: solution.npy")
    if not (solution.get("file") or solution.get("values") is not None):
        raise ConfigError(
            "solution needs a converged solution vector",
            "solution:\n  file: solution.npy    # or:  values: [0.0, 1.0, ...]")

    unknowns = _resolve_unknowns(declared, solution, base_dir)

    sens = _require(raw, "sensitivity", "sensitivity:\n  order: 1\n  backend: otilib")
    order = _require(sens, "order", "sensitivity:\n  order: 1", where="sensitivity")
    try:
        order = int(order)
    except (TypeError, ValueError):
        raise ConfigError("sensitivity.order must be an integer >= 1",
                          "sensitivity:\n  order: 1")
    if order < 1:
        raise ConfigError("sensitivity.order must be >= 1 (got %d)" % order,
                          "sensitivity:\n  order: 1")
    backend = str(sens.get("backend", "otilib")).lower()

    tangent = raw.get("tangent") or {}
    if rtype == "executable" and not tangent:
        # black-box may return the tangent in its response; a file is also fine
        tangent = {"type": "response"}

    cfg = UserConfig(
        raw=raw, base_dir=base_dir, name=str(name), unknowns=unknowns,
        residual_type=rtype, parameters=parameters, order=order, backend=backend,
        residual=residual, tangent=tangent, solution=solution,
        constraints=raw.get("constraints") or {}, state=raw.get("state") or {},
        time=raw.get("time") or {},
        output_dir=(raw.get("output") or {}).get("dir", "resasm_output"),
        validation=raw.get("validation") or {})
    return cfg


def solution_size(solution: Dict[str, Any], base_dir: str) -> Optional[int]:
    """Length of the solution vector, or None if it cannot be determined yet
    (e.g. the .npy is not on disk). Never raises."""
    vals = solution.get("values")
    if vals is not None:
        try:
            return int(np.asarray(vals, float).ravel().size)
        except Exception:                           # noqa: BLE001
            return None
    f = solution.get("file")
    if not f:
        return None
    p = f if os.path.isabs(f) else os.path.join(base_dir, f)
    if not os.path.exists(p):
        return None
    try:
        return int(np.asarray(np.load(p, allow_pickle=True), float).ravel().size)
    except Exception:                               # noqa: BLE001
        return None


def _resolve_unknowns(declared: Optional[int], solution: Dict[str, Any],
                      base_dir: str) -> int:
    """problem.unknowns is optional -> infer from the solution vector.

    - declared only            -> use it
    - solution only            -> infer (this is the low-burden path)
    - both                     -> must agree, else a clear ConfigError
    - neither determinable     -> ask for one of them
    """
    found = solution_size(solution, base_dir)

    if declared is not None and found is not None:
        if declared != found:
            raise ConfigError(
                "problem.unknowns is %d but the solution vector has %d entries"
                % (declared, found),
                "problem:\n  unknowns: %d      # or delete this line and let it be "
                "inferred" % found)
        return declared
    if declared is not None:
        return declared
    if found is not None:
        return found
    raise ConfigError(
        "cannot determine the number of unknowns: the solution vector is not "
        "readable yet and problem.unknowns is not set",
        "problem:\n  unknowns: 1        # or make solution.file readable\n"
        "solution:\n  file: solution.npy")


def _validate_residual(residual: Dict[str, Any], rtype: str) -> None:
    if rtype == "python":
        if not residual.get("module"):
            raise ConfigError(
                "residual.module is required for a Python residual",
                "residual:\n  type: python\n  module: user_residual.py\n"
                "  function: residual")
        if not residual.get("function"):
            residual["function"] = "residual"
    elif rtype == "executable":
        if not residual.get("command"):
            raise ConfigError(
                "residual.command is required for a black-box executable",
                "residual:\n  type: executable\n"
                "  command: ./my_solver --request {request} --response {response}")
    elif rtype == "element":
        if not residual.get("module"):
            raise ConfigError(
                "residual.module is required for an element residual provider",
                "residual:\n  type: element\n  module: user_elements.py")


def _require(d: Any, key: str, fix: str, where: Optional[str] = None) -> Any:
    if not isinstance(d, dict) or key not in d or d[key] is None:
        loc = "%s.%s" % (where, key) if where else key
        raise ConfigError("Missing required field: %s" % loc, fix)
    return d[key]
