"""`resasm init` — an interactive wizard that writes a minimal resasm.yml.

Non-interactive callers can pass an ``answers`` dict (used by tests).
"""

from __future__ import annotations

import os
from typing import Any, Callable, Dict, Optional


_PATHS = {"1": "python", "2": "cpp", "3": "blackbox",
          "python": "python", "cpp": "cpp", "compiled": "cpp",
          "blackbox": "blackbox", "executable": "blackbox"}


def init_wizard(out_path: str = "resasm.yml",
                answers: Optional[Dict[str, Any]] = None,
                prompt: Callable[[str], str] = input,
                echo: Callable[[str], None] = print) -> str:
    """Collect a few answers and write resasm.yml. Returns the written text."""
    a = answers or {}

    def ask(key, question, default=None):
        if key in a:
            return a[key]
        suffix = " [%s]" % default if default is not None else ""
        val = prompt("%s%s: " % (question, suffix)).strip()
        return val or default

    path = _PATHS.get(str(ask("path", "Residual path: 1) Python  2) C++/Fortran  "
                                       "3) black-box executable", "1")).lower(), "python")
    name = ask("name", "Problem name", "my_job")
    unknowns = int(ask("unknowns", "Number of unknowns (DOFs)", 1))
    raw_params = ask("parameters", "Parameters as name=value,comma-separated",
                     "k=2.0,f=16.0")
    parameters = _parse_params(raw_params)
    solution = ask("solution", "Solution vector file", "solution.npy")
    order = int(ask("order", "Derivative order", 1))
    backend = ask("backend", "Backend (otilib)", "otilib")

    if path == "python":
        residual = {"type": "python", "module": ask("module", "Residual module",
                                                     "user_residual.py"),
                    "function": "residual"}
        tangent = {"type": "python", "function": "tangent"}
    elif path == "blackbox":
        residual = {"type": "executable",
                    "command": ask("command", "Executable command",
                                   "python my_solver.py --request {request} "
                                   "--response {response}")}
        tangent = {"type": ask("tangent_type", "Tangent source (file/response)",
                               "response")}
        if tangent["type"] == "file":
            tangent["file"] = ask("tangent_file", "Tangent file", "tangent.npz")
    else:  # cpp/fortran compiled -> exposed as an executable wrapper
        residual = {"type": "executable",
                    "command": ask("command", "Executable command",
                                   "./my_solver --request {request} "
                                   "--response {response}")}
        tangent = {"type": "response"}

    text = _emit_yaml(name, unknowns, residual, tangent, parameters, solution,
                      order, backend)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(text)
    echo("wrote %s" % out_path)
    return text


def _parse_params(raw) -> Dict[str, float]:
    if isinstance(raw, dict):
        return {k: float(v) for k, v in raw.items()}
    out = {}
    for item in str(raw).split(","):
        if "=" in item:
            k, v = item.split("=", 1)
            try:
                out[k.strip()] = float(v)
            except ValueError:
                out[k.strip()] = 0.0
    return out or {"k": 1.0}


def _emit_yaml(name, unknowns, residual, tangent, parameters, solution, order,
               backend) -> str:
    lines = ["problem:", "  name: %s" % name, "  unknowns: %d" % unknowns, "",
             "residual:", "  type: %s" % residual["type"]]
    for k in ("module", "function", "command"):
        if residual.get(k):
            lines.append("  %s: %s" % (k, residual[k]))
    lines += ["", "tangent:", "  type: %s" % tangent["type"]]
    for k in ("function", "file"):
        if tangent.get(k):
            lines.append("  %s: %s" % (k, tangent[k]))
    lines += ["", "parameters:"]
    for k, v in parameters.items():
        lines.append("  %s: %s" % (k, v))
    lines += ["", "solution:", "  file: %s" % solution, "",
              "sensitivity:", "  order: %d" % order, "  backend: %s" % backend, ""]
    return "\n".join(lines)
