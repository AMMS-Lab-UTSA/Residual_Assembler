"""Run an assembly recipe (Path A): build R from ingredients, then sensitivities.

    ingredients -> Residual_Assembler assembles R -> OTI overloads parameters
                -> extract R^(p) -> solve T U^(p) = -R^(p) -> package

This module only *orchestrates*. The assembly and the OTI algebra live in
`residual_core` (unchanged); the output package is written by the same
`resasm_user.output.write_outputs` used by every other path, so an assembly job
produces exactly the same private/ + public/ deliverable.

Honesty gate
------------
Assembling R and DIFFERENTIATING R are different capabilities. A backend whose
kernel uses numpy float arrays can do the first and not the second (it cannot
carry an OTI number). We refuse to run a sensitivity through such a backend and
say precisely why, rather than emit a plausible-looking wrong number.
"""

from __future__ import annotations

import os
import time as _time
from typing import Any, Dict

import numpy as np

from .config import ConfigError
from .recipe import Recipe, load_recipe, describe, sensitivity_capability


def _problem(r: Recipe):
    """Build the ResidualProblem from the recipe's mesh."""
    from residual_core import ResidualProblem
    path = r.path(r.mesh.file)
    if r.mesh.format == "neutral_json":
        return ResidualProblem.from_neutral(path)
    if r.mesh.format == "abaqus_inp":
        return ResidualProblem.from_abaqus(path)
    raise ConfigError("unknown mesh format %r (use a .inp or .json mesh)"
                      % r.mesh.format)


def check_recipe(path: str) -> "CheckResult":
    """`resasm check` for an assembly recipe: report given / inferred / missing."""
    r = load_recipe(path)
    text = describe(r)
    ok = not r.missing and not r.blockers
    return CheckResult(ok=ok, recipe=r, text=text)


class CheckResult(object):
    def __init__(self, ok, recipe, text):
        self.ok = ok
        self.recipe = recipe
        self.text = text

    def render(self):
        return self.text


def run_recipe(path: str):
    """`resasm run` for an assembly recipe."""
    from .runner import RunResult
    from .output import write_outputs

    t0 = _time.time()
    r = load_recipe(path)

    if r.missing or r.blockers:
        raise ConfigError(
            "the assembly recipe is not complete:\n\n" + describe(r))

    cap = sensitivity_capability(r)
    if not cap["oti_differentiable"]:
        raise ConfigError(
            "cannot compute sensitivities for this model.\n\n"
            "  assemble R          : yes\n"
            "  OTI-differentiate R : NO\n\n"
            "%s\n\n"
            "What you CAN do today:\n"
            "  * assemble and verify the residual  (resasm assemble / resasm verify)\n"
            "  * use the black-box path: your own solver returns the residual\n"
            "    coefficients      (resasm init --template blackbox-order2)\n"
            "  * plug in an OTI-transformed UMAT (companion project) once available"
            % cap.get("reason", ""))

    prob = _problem(r)

    # solution field (converged u)
    upath = r.path(r.fields.solution)
    if upath and os.path.exists(upath):
        prob.set_solution(np.load(upath, allow_pickle=True).ravel())

    order = int(r.sensitivity.get("order", 1))
    # assembly sub-mode: formulation (self-contained elements) | material-replay |
    # stress-driven. Default 'formulation'.
    mode = str(r.sensitivity.get("mode") or "formulation")

    pkg = prob.sensitivity_package(mode=mode, parameters=list(r.parameters),
                                   max_order=order, algebra="OTI",
                                   generate_rhs=True, backend="otilib")
    if not getattr(pkg, "runnable", False):
        raise ConfigError("cannot assemble in mode '%s': %s"
                          % (mode, getattr(pkg, "minimum_missing", "unknown")))

    sens = pkg.sensitivity
    res = pkg.residual
    names = [n for n, _ in sorted(sens.parameter_map.items(), key=lambda kv: kv[1])]
    ndof = int(sens.ndof)
    free = res.free_mask if res.free_mask is not None else np.ones(ndof, bool)
    T = res.tangent.T

    from . import oti_global
    labels = oti_global.direction_labels(len(names), order)

    result = {
        "parameter_names": names, "ndof": ndof, "free_mask": np.asarray(free, bool),
        "R_real": np.asarray(res.R, float),
        "tangent": np.asarray(T, float), "tangent_source": res.tangent.source,
        "R_orders": {p: np.asarray(a, float) for p, a in sens._R.items()},
        "U_orders": {p: np.asarray(a, float) for p, a in sens._solved_U.items()},
        "direction_labels": labels,
        "validation": _validation(res, pkg),
        "diagnostics": {"engine": "otilib", "mode": mode,
                        "assembly": res.diagnostics},
        "timing": {"total_seconds": round(_time.time() - t0, 4),
                   "orders": order, "unknowns": ndof},
    }

    class _Cfg(object):
        name = r.name
        order = int(r.sensitivity.get("order", 1))
        backend = str(r.sensitivity.get("backend", "otilib"))
        residual_type = "assemble"

    out_dir = r.path(r.output_dir)
    paths = write_outputs(out_dir, _Cfg(), result)
    return RunResult(ok=True, output_dir=out_dir, private_dir=paths["private"],
                     public_dir=paths["public"],
                     summary={"mode": "assemble", "parameters": names,
                              "order": order,
                              "tangent_source": res.tangent.source,
                              "residual_free_norm":
                                  result["validation"].get("residual_free_norm"),
                              "orders_solved": sorted(result["U_orders"])})


def _validation(res, pkg) -> Dict[str, Any]:
    v: Dict[str, Any] = {"status": "ok"}
    R = np.asarray(res.R, float)
    free = res.free_mask if res.free_mask is not None else np.ones(R.size, bool)
    v["residual_free_norm"] = float(np.linalg.norm(R[np.asarray(free, bool)]))
    rep = getattr(pkg, "validation", None)
    if rep is not None and getattr(rep, "checks", None):
        v["assembly_checks"] = [
            {"name": c.name, "quantity": c.quantity, "value": c.value,
             "passed": c.passed, "detail": c.detail} for c in rep.checks]
    return v
