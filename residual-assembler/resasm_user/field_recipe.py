"""Field-driven residual-sensitivity recipe (Milestone M2).

Makes the M1 engine (``residual_core.core.field_sensitivity``) runnable
end-to-end from user files::

    resasm run sensitivity.yaml

where ``sensitivity.yaml`` is::

    analysis:
      type: field_residual_sensitivity
      kinematics: small_strain
    mesh:
      format: abaqus            # or 'neutral'
      file: model.inp
    results:
      format: resasm_derivative_fields_v1
      file: derivative_fields.json
    parameters: [E, nu]
    assumptions:
      parameter_independent_geometry: true
      parameter_independent_loads: true
      parameter_independent_boundaries: true
    output:
      directory: results
      save_tangent: true
      save_residual_derivatives: true
      save_displacement_sensitivities: true

The runner loads the model, loads + validates the derivative fields (whose
metadata is the single source of truth for the SDV layout), splits them into
per-IP tangent and stress-derivative fields, solves K U_,a = -R_,a for every
parameter at once, and writes ``results/{displacement_sensitivities.npz,
residual_derivatives.npz, tangent.npz, manifest.json}``.

It changes NONE of the M1 mathematics; it only connects inputs to the engine and
guards the supported scope (small-strain C3D8, parameter-independent geometry /
loads / BCs), refusing anything else with a clear message instead of a partial
answer.
"""

from __future__ import annotations

import os
import time as _time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from .config import ConfigError
from . import _miniyaml as _yaml

ANALYSIS_TYPE = "field_residual_sensitivity"
MANIFEST_SCHEMA = "resasm_field_sensitivity_manifest_v1"
_DEFAULT_SAVE = {"save_displacement_sensitivities": True,
                 "save_residual_derivatives": True,
                 "save_tangent": True}


# --------------------------------------------------------------------------- #
# recipe object
# --------------------------------------------------------------------------- #
@dataclass
class FieldRecipe:
    name: str = "field_sensitivity"
    kinematics: str = "small_strain"
    mesh_file: Optional[str] = None
    mesh_format: Optional[str] = None            # abaqus | neutral
    results_file: Optional[str] = None
    results_format: Optional[str] = None
    parameters: List[str] = field(default_factory=list)
    assumptions: Dict[str, Any] = field(default_factory=dict)
    output_dir: str = "results"
    save: Dict[str, bool] = field(default_factory=lambda: dict(_DEFAULT_SAVE))
    base_dir: str = "."

    def path(self, rel: Optional[str]) -> Optional[str]:
        if not rel:
            return None
        return rel if os.path.isabs(rel) else os.path.join(self.base_dir, rel)


# --------------------------------------------------------------------------- #
# detection + loading
# --------------------------------------------------------------------------- #
def _analysis_type_of(raw: Any) -> Optional[str]:
    if not isinstance(raw, dict):
        return None
    analysis = raw.get("analysis")
    if isinstance(analysis, dict):
        return str(analysis.get("type") or "").strip() or None
    if isinstance(analysis, str):
        return analysis.strip() or None
    return None


def is_field_recipe(path: str) -> bool:
    """True if the config at ``path`` declares analysis.type = field_residual_sensitivity."""
    try:
        raw = _yaml.load_file(path)
    except Exception:                                    # noqa: BLE001
        return False
    return _analysis_type_of(raw) == ANALYSIS_TYPE


def _mesh_format(fmt: Optional[str], file: Optional[str]) -> Optional[str]:
    if fmt:
        f = str(fmt).lower()
        if f in ("abaqus", "abaqus_inp", "inp"):
            return "abaqus"
        if f in ("neutral", "neutral_json", "json"):
            return "neutral"
        return f
    if file:
        ext = os.path.splitext(file)[1].lower()
        return {".inp": "abaqus", ".json": "neutral"}.get(ext)
    return None


def load_field_recipe(path: str) -> FieldRecipe:
    """Parse and validate a field_residual_sensitivity recipe."""
    if not os.path.exists(path):
        raise ConfigError("recipe file not found: %s" % path)
    raw = _yaml.load_file(path)
    if not isinstance(raw, dict):
        raise ConfigError("recipe must be a mapping (key: value pairs)")
    if _analysis_type_of(raw) != ANALYSIS_TYPE:
        raise ConfigError("analysis.type must be %r for the field-driven path, got %r"
                          % (ANALYSIS_TYPE, _analysis_type_of(raw)))

    r = FieldRecipe()
    r.base_dir = os.path.dirname(os.path.abspath(path))

    analysis = raw.get("analysis") or {}
    if isinstance(analysis, dict):
        r.kinematics = str(analysis.get("kinematics") or "small_strain").lower()
    r.name = str(raw.get("name") or (analysis.get("name") if isinstance(analysis, dict) else None)
                 or "field_sensitivity")

    mesh = raw.get("mesh") or {}
    if isinstance(mesh, str):
        mesh = {"file": mesh}
    r.mesh_file = mesh.get("file")
    r.mesh_format = _mesh_format(mesh.get("format"), r.mesh_file)

    results = raw.get("results") or {}
    if isinstance(results, str):
        results = {"file": results}
    r.results_file = results.get("file")
    r.results_format = results.get("format")

    params = raw.get("parameters")
    if isinstance(params, dict):
        r.parameters = [str(k) for k in params]
    elif isinstance(params, list):
        r.parameters = [str(p) for p in params]
    elif isinstance(params, str):
        r.parameters = [params]

    assumptions = raw.get("assumptions") or {}
    r.assumptions = dict(assumptions) if isinstance(assumptions, dict) else {}

    out = raw.get("output") or {}
    if isinstance(out, str):
        out = {"directory": out}
    r.output_dir = out.get("directory") or out.get("dir") or "results"
    for k in _DEFAULT_SAVE:
        if k in out:
            r.save[k] = bool(out[k])

    _validate_recipe(r)
    return r


def _validate_recipe(r: FieldRecipe) -> None:
    if not r.mesh_file:
        raise ConfigError("recipe is missing mesh.file")
    if r.mesh_format not in ("abaqus", "neutral"):
        raise ConfigError("mesh.format must be 'abaqus' or 'neutral' (got %r)"
                          % r.mesh_format)
    if not r.results_file:
        raise ConfigError("recipe is missing results.file (the derivative-field JSON)")
    if r.results_format and r.results_format != "resasm_derivative_fields_v1":
        raise ConfigError("results.format must be 'resasm_derivative_fields_v1', got %r"
                          % r.results_format)
    if r.kinematics != "small_strain":
        raise ConfigError(
            "kinematics %r is not supported; the field-driven path is small_strain "
            "only (finite strain is a later milestone)" % r.kinematics)
    if not r.parameters:
        raise ConfigError("recipe is missing 'parameters' (the ordered list of "
                          "parameter names to differentiate, e.g. [E, nu])")


# --------------------------------------------------------------------------- #
# scope guards
# --------------------------------------------------------------------------- #
_PARAM_DEP_ASSUMPTIONS = ("parameter_independent_geometry",
                          "parameter_independent_loads",
                          "parameter_independent_boundaries")


def _guard_assumptions(r: FieldRecipe) -> None:
    for key in _PARAM_DEP_ASSUMPTIONS:
        if key in r.assumptions and not bool(r.assumptions[key]):
            what = key.replace("parameter_independent_", "")
            raise ConfigError(
                "assumptions.%s is false: parameter-dependent %s is not supported "
                "by the field-driven path (the residual method here assumes "
                "d f_ext/d a = 0 and d u_prescribed/d a = 0)." % (key, what))


def _guard_model(r: FieldRecipe, model, abaqus_model) -> None:
    # unsupported element types
    bad = sorted({e.etype.upper() for e in model.elements.values()
                  if e.etype.upper() != "C3D8"})
    if bad:
        raise ConfigError("unsupported element type(s) %s; the field-driven path "
                          "supports C3D8 only" % ", ".join(bad))
    if not any(e.etype.upper() == "C3D8" for e in model.elements.values()):
        raise ConfigError("the mesh has no C3D8 elements to assemble")

    # *Equation constraints (parsed but never applied -> refuse rather than ignore)
    if getattr(model, "equations", None):
        raise ConfigError("*Equation (linear MPC) constraints are present but not "
                          "supported by the field-driven path; remove them or use a "
                          "model without linear constraints")

    if abaqus_model is not None:
        if getattr(abaqus_model, "dsloads", None):
            raise ConfigError("*Dsload / distributed surface loads are present but "
                              "not supported by the field-driven path")
        unsup = [k for k in getattr(abaqus_model, "unsupported_keywords", []) or []
                 if any(t in k.lower() for t in ("contact", "interaction", "dload",
                                                 "dflux", "tie"))]
        if unsup:
            raise ConfigError("unsupported model feature(s): %s (contact / distributed "
                              "load / tie are out of scope for the field-driven path)"
                              % ", ".join(unsup))


# --------------------------------------------------------------------------- #
# the runner
# --------------------------------------------------------------------------- #
def run_field_recipe(path: str):
    """Execute a field_residual_sensitivity recipe end-to-end.

    Returns a :class:`resasm_user.runner.RunResult`.
    """
    from .runner import RunResult
    from residual_core import ResidualProblem
    from residual_core.core.dof_manager import DofManager
    from residual_core.core.field_sensitivity import (
        solve_field_sensitivities, fields_from_statev, FieldSensitivityError)
    from residual_core.io.derivative_fields import (
        DerivativeFieldError, displacements_to_vector)

    t0 = _time.time()
    r = load_field_recipe(path)
    _guard_assumptions(r)

    mesh_path = r.path(r.mesh_file)
    if not os.path.exists(mesh_path):
        raise ConfigError("mesh file not found: %s" % mesh_path)
    results_path = r.path(r.results_file)
    if not os.path.exists(results_path):
        raise ConfigError("derivative-field file not found: %s" % results_path)

    # 1) load model (exercises the wired ResidualProblem.attach_results path)
    if r.mesh_format == "neutral":
        prob = ResidualProblem.from_neutral(mesh_path)
    else:
        prob = ResidualProblem.from_abaqus(mesh_path)
    try:
        prob.attach_results(results_path)
    except DerivativeFieldError as exc:
        raise ConfigError(str(exc))
    if prob.results is None:
        raise ConfigError("results file %s is not a resasm_derivative_fields_v1 export"
                          % results_path)
    df = prob.results
    model = prob.model

    # 2) scope guards
    _guard_model(r, model, getattr(prob, "_abaqus", None))
    if r.kinematics != df.kinematics:
        raise ConfigError("recipe kinematics %r disagrees with the field file's %r"
                          % (r.kinematics, df.kinematics))
    layout_params = df.parameters
    unknown = [p for p in r.parameters if p not in layout_params]
    if unknown:
        raise ConfigError(
            "parameter(s) %r are not in the derivative-field metadata (available: "
            "%r); a misspelled parameter must never become a silent zero sensitivity"
            % (unknown, layout_params))

    # 3) build the global solution vector + split the fields, then solve
    dm = DofManager(model.nodes.keys())
    u_vec = displacements_to_vector(df.displacements, dm)
    try:
        tangent_fields, stress_fields = fields_from_statev(
            df.statev, df.sdv_layout, parameters=r.parameters)
        res = solve_field_sensitivities(
            model=model, solution=u_vec, tangent_fields=tangent_fields,
            stress_derivative_fields=stress_fields, parameters=r.parameters,
            dof_manager=dm, integration=df.integration)
    except FieldSensitivityError as exc:
        raise ConfigError(str(exc))

    # 4) write outputs
    out_dir = r.path(r.output_dir)
    paths = _write_outputs(out_dir, res, r, df, mesh_path, results_path,
                           round(_time.time() - t0, 4))

    summary = {
        "analysis": ANALYSIS_TYPE,
        "parameters": res.parameters,
        "ndof": res.ndof,
        "n_elements": res.diagnostics.get("n_elements"),
        "tangent_source": res.diagnostics.get("tangent_source"),
        "outputs": paths["files"],
        "du_da_norms": {p: float(np.linalg.norm(res.du_da(p))) for p in res.parameters},
    }
    return RunResult(ok=True, output_dir=out_dir, private_dir=out_dir,
                     public_dir=out_dir, summary=summary)


# --------------------------------------------------------------------------- #
# output
# --------------------------------------------------------------------------- #
def _write_outputs(out_dir, res, recipe, df, mesh_path, results_path, seconds):
    import json
    os.makedirs(out_dir, exist_ok=True)
    files: List[str] = []

    if recipe.save.get("save_displacement_sensitivities", True):
        named = {"du_d%s" % p: res.du_da(p) for p in res.parameters}
        np.savez(os.path.join(out_dir, "displacement_sensitivities.npz"),
                 U=res.displacement_sensitivities,
                 parameters=np.array(res.parameters, dtype=object),
                 free_mask=res.free_mask, prescribed_idx=res.prescribed_idx, **named)
        files.append("displacement_sensitivities.npz")

    if recipe.save.get("save_residual_derivatives", True):
        named = {"dR_d%s" % p: res.R_da(p) for p in res.parameters}
        np.savez(os.path.join(out_dir, "residual_derivatives.npz"),
                 R=res.residual_derivatives,
                 parameters=np.array(res.parameters, dtype=object), **named)
        files.append("residual_derivatives.npz")

    if recipe.save.get("save_tangent", True):
        np.savez(os.path.join(out_dir, "tangent.npz"), K=res.K)
        files.append("tangent.npz")

    manifest = {
        "schema": MANIFEST_SCHEMA,
        "analysis": ANALYSIS_TYPE,
        "created": _time.strftime("%Y-%m-%dT%H:%M:%S"),
        "runtime_seconds": seconds,
        "parameters": res.parameters,
        "ndof": res.ndof,
        "n_free": res.diagnostics.get("n_free"),
        "n_prescribed": res.diagnostics.get("n_prescribed"),
        "n_elements": res.diagnostics.get("n_elements"),
        "tangent_source": res.diagnostics.get("tangent_source"),
        "kinematics": df.kinematics,
        "integration": res.diagnostics.get("integration"),
        "element_type": df.element_type,
        "voigt_order": df.metadata.get("voigt_order"),
        "integration_point_order": df.metadata.get("integration_point_order"),
        "sdv_indexing": df.metadata.get("sdv_indexing"),
        "sdv_layout": df.sdv_layout,
        "assumptions": recipe.assumptions,
        "inputs": {"mesh": os.path.basename(mesh_path),
                   "derivative_fields": os.path.basename(results_path)},
        "outputs": files,
        "du_da_norms": {p: float(np.linalg.norm(res.du_da(p))) for p in res.parameters},
        "prescribed_dof_sensitivity_is_zero": bool(
            np.all(res.displacement_sensitivities[res.prescribed_idx] == 0.0)),
    }
    with open(os.path.join(out_dir, "manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    files.append("manifest.json")

    _write_summary_md(os.path.join(out_dir, "summary.md"), res, recipe, df, manifest)
    files.append("summary.md")
    return {"dir": out_dir, "files": files}


def _write_summary_md(path, res, recipe, df, manifest):
    lines = [
        "# Field-driven residual sensitivity — %s" % recipe.name,
        "",
        "| item | value |", "|---|---|",
        "| analysis | %s |" % ANALYSIS_TYPE,
        "| element type / kinematics | %s / %s |" % (df.element_type, df.kinematics),
        "| unknowns (DOFs) | %d |" % res.ndof,
        "| free / prescribed | %d / %d |"
        % (manifest["n_free"], manifest["n_prescribed"]),
        "| elements | %d |" % manifest["n_elements"],
        "| tangent source | %s |" % manifest["tangent_source"],
        "| parameters | %s |" % ", ".join(res.parameters),
        "",
        "## du/da norms", "",
        "| parameter | ‖du/da‖ |", "|---|---|",
    ]
    for p in res.parameters:
        lines.append("| %s | %.6e |" % (p, float(np.linalg.norm(res.du_da(p)))))
    lines += ["",
              "> Prescribed-DOF sensitivities are exactly zero: %s."
              % manifest["prescribed_dof_sensitivity_is_zero"],
              "> Files: displacement_sensitivities.npz, residual_derivatives.npz, "
              "tangent.npz, manifest.json.", ""]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
