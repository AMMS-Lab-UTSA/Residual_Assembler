"""ResidualAssemblyRecipe — the ingredients we need to BUILD your residual.

The central object of the primary (Path A) workflow.

    The user should not provide R.
    The user should provide enough ingredients for us to build R.

Abaqus hides the global residual. But R is not magic:

    R(u, a) = F_internal(u, a, q) - F_external(a, t) + F_constraints(u, t)
    R_e     = INT_Omega_e  B^T sigma(u, a, q) dOmega  -  f_e^ext

so we ask for the things that go INTO it — mesh, element formulation, material,
solution field, state/history, loads/BCs — not for R itself.

This module does three jobs:

  1. **Model the recipe** (`ResidualAssemblyRecipe`) — the analog of the old UMAT
     project's transformation contract.
  2. **Progressive inference** — infer everything inferable (element type, ndof,
     integration rule, DOF map, backends, BCs) so the user types as little as
     possible.
  3. **Honest capability gating** — say exactly what is missing, and refuse to
     pretend a sensitivity is available when the backend cannot carry an OTI
     number (see `Formulation.oti_differentiable`).

It performs NO OTI algebra and NO assembly itself: it resolves ingredients and
then hands off to `residual_core` (`ResidualProblem`).
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from .config import ConfigError
from . import _miniyaml as _yaml


# --------------------------------------------------------------------------- #
# Ingredient objects
# --------------------------------------------------------------------------- #
@dataclass
class MeshSpec:
    file: Optional[str] = None
    format: Optional[str] = None        # abaqus_inp | neutral_json
    n_nodes: Optional[int] = None
    n_elements: Optional[int] = None
    element_types: List[str] = field(default_factory=list)


@dataclass
class FormulationSpec:
    """Element formulation per element type. Normally INFERRED from the mesh."""
    backends: Dict[str, str] = field(default_factory=dict)   # etype -> backend key
    element: Optional[str] = None                            # user override
    kinematics: Optional[str] = None
    integration: Optional[str] = None                        # implied by the backend


@dataclass
class MaterialSpec:
    type: Optional[str] = None          # umat | elastic | builtin | none
    file: Optional[str] = None          # UMAT source
    name: Optional[str] = None
    properties: Dict[str, float] = field(default_factory=dict)


@dataclass
class FieldsSpec:
    solution: Optional[str] = None      # converged U (.npy)
    element_fields: Optional[str] = None  # exported IP stress etc. (fields.json)
    source: Optional[str] = None        # e.g. odb_export


@dataclass
class StimuliSpec:
    loads: Optional[str] = None
    boundary_conditions: Optional[str] = None
    from_mesh: bool = True              # BCs/loads usually live in the .inp


@dataclass
class Recipe:
    """The full residual-assembly recipe."""
    name: str = "assembly_job"
    mode: str = "assemble"
    mesh: MeshSpec = field(default_factory=MeshSpec)
    formulation: FormulationSpec = field(default_factory=FormulationSpec)
    material: MaterialSpec = field(default_factory=MaterialSpec)
    fields: FieldsSpec = field(default_factory=FieldsSpec)
    stimuli: StimuliSpec = field(default_factory=StimuliSpec)
    state: Dict[str, Any] = field(default_factory=dict)
    time: Dict[str, Any] = field(default_factory=dict)
    constraints: Dict[str, Any] = field(default_factory=dict)
    parameters: List[str] = field(default_factory=list)      # ["mat.key", ...]
    tangent: Dict[str, Any] = field(default_factory=dict)
    sensitivity: Dict[str, Any] = field(default_factory=dict)
    output_dir: str = "resasm_output"
    base_dir: str = "."

    # ---- filled in by inference ------------------------------------------
    inferred: Dict[str, str] = field(default_factory=dict)   # what -> where from
    missing: List[str] = field(default_factory=list)
    blockers: List[str] = field(default_factory=list)        # hard "cannot do it"

    def path(self, rel: Optional[str]) -> Optional[str]:
        if not rel:
            return None
        return rel if os.path.isabs(rel) else os.path.join(self.base_dir, rel)


# --------------------------------------------------------------------------- #
# Loading + progressive inference
# --------------------------------------------------------------------------- #
def _as_dict(v, key="file"):
    """Accept both `mesh: model.inp` and `mesh: {file: model.inp}`."""
    if v is None:
        return {}
    if isinstance(v, dict):
        return v
    return {key: v}


def is_assembly_config(raw: Dict[str, Any]) -> bool:
    """An assembly recipe is anything that names a mesh (and does not hand us a
    ready-made residual)."""
    if not isinstance(raw, dict):
        return False
    if str(raw.get("mode", "")).lower() in ("assemble", "assembly"):
        return True
    return "mesh" in raw and "residual" not in raw


def load_recipe(path: str) -> Recipe:
    """Parse an assembly `resasm.yml` and infer everything inferable."""
    if not os.path.exists(path):
        raise ConfigError("Config file not found: %s" % path)
    raw = _yaml.load_file(path)
    if not isinstance(raw, dict):
        raise ConfigError("Config must be a mapping (key: value pairs)")

    r = Recipe()
    r.base_dir = os.path.dirname(os.path.abspath(path))
    prob = _as_dict(raw.get("problem"), "name")
    r.name = str(prob.get("name") or raw.get("name") or "assembly_job")
    r.mode = str(raw.get("mode") or "assemble").lower()

    # mesh: `mesh: model.inp`  or  `mesh: {file: ..., format: ...}`
    mesh = _as_dict(raw.get("mesh"))
    r.mesh.file = mesh.get("file")
    r.mesh.format = mesh.get("format")

    # fields: `solution: U.npy`  or  `fields: {solution: {file: U.npy}}`
    flds = raw.get("fields") or {}
    sol = flds.get("solution") if isinstance(flds, dict) else None
    r.fields.solution = (_as_dict(sol).get("file") if sol is not None
                         else (raw.get("solution") if isinstance(raw.get("solution"), str)
                               else _as_dict(raw.get("solution")).get("file")))
    ef = flds.get("element_fields") if isinstance(flds, dict) else None
    r.fields.element_fields = (_as_dict(ef).get("file") if ef is not None
                               else raw.get("element_fields"))
    if isinstance(sol, dict):
        r.fields.source = sol.get("source")

    # material: `material: umat.f`  or  a full block
    mat = _as_dict(raw.get("material"), "file")
    r.material.file = mat.get("file")
    r.material.type = mat.get("type")
    r.material.name = mat.get("name")
    r.material.properties = dict(mat.get("properties") or {})

    form = _as_dict(raw.get("formulation"), "element")
    r.formulation.element = form.get("element")
    r.formulation.kinematics = form.get("kinematics")
    r.formulation.integration = form.get("integration")

    stim = raw.get("stimuli") or {}
    if isinstance(stim, dict):
        r.stimuli.loads = _as_dict(stim.get("loads")).get("file")
        r.stimuli.boundary_conditions = _as_dict(
            stim.get("boundary_conditions")).get("file")

    r.state = raw.get("state") or {}
    r.time = raw.get("time") or {}
    r.constraints = raw.get("constraints") or {}
    r.tangent = raw.get("tangent") or {}
    r.sensitivity = raw.get("sensitivity") or {}
    r.output_dir = (raw.get("output") or {}).get("dir", "resasm_output")

    params = raw.get("parameters")
    if isinstance(params, dict):
        # {"spring.k": 2.0} or {"E": 210000.0} -> names (values live in the model)
        r.parameters = [str(k) for k in params]
    elif isinstance(params, list):
        r.parameters = [str(p) for p in params]

    infer(r)
    return r


def infer(r: Recipe) -> Recipe:
    """Fill in everything derivable, and record WHERE it came from.

    Inferred: mesh format, element types, per-type backend, node/element counts,
    integration rule, sensitivity order/backend, output dir.
    Never guessed: the solution field, the material evaluator, the parameters.
    """
    # --- mesh format from the extension
    if r.mesh.file and not r.mesh.format:
        ext = os.path.splitext(r.mesh.file)[1].lower()
        r.mesh.format = {".inp": "abaqus_inp", ".json": "neutral_json"}.get(ext)
        if r.mesh.format:
            r.inferred["mesh.format"] = "from the file extension '%s'" % ext

    if not r.mesh.file:
        r.missing.append("mesh: the model/mesh file (e.g. model.inp)")
        return r

    # --- parse the mesh: element types, counts, backends
    try:
        model = _load_model(r)
    except Exception as exc:                              # noqa: BLE001
        r.missing.append("mesh could not be read: %s" % exc)
        return r

    etypes = {}
    for e in model.elements.values():
        etypes[e.etype.upper()] = etypes.get(e.etype.upper(), 0) + 1
    r.mesh.element_types = sorted(etypes)
    r.mesh.n_nodes = len(model.nodes)
    r.mesh.n_elements = len(model.elements)
    r.inferred["mesh.element_types"] = "from %s (%s)" % (
        os.path.basename(r.mesh.file),
        ", ".join("%s x%d" % (k, v) for k, v in sorted(etypes.items())))
    r.inferred["dof_map"] = "built from the mesh (%d nodes)" % r.mesh.n_nodes

    # --- backend per element type (the element formulation)
    from residual_core.formulations import default_formulations
    reg = default_formulations()
    bound = dict(getattr(model, "element_formulation", {}) or {})
    for et in r.mesh.element_types:
        picked = r.formulation.element
        if not picked:
            # honour an explicit binding in the model, else the default policy
            for eid, el in model.elements.items():
                if el.etype.upper() == et and eid in bound:
                    picked = bound[eid]
                    break
        if not picked:
            from residual_core.core.model import default_formulation_policy
            picked = default_formulation_policy(et)
        if picked and picked in reg:
            r.formulation.backends[et] = picked
        else:
            r.blockers.append(
                "no element formulation for '%s' — not implemented "
                "(supported: %s)" % (et, ", ".join(sorted(
                    t for f in reg.values() for t in f.element_types))))
    if r.formulation.backends:
        r.inferred["formulation.backend"] = "auto-selected per element type: " + \
            ", ".join("%s -> %s" % (k, v) for k, v in r.formulation.backends.items())
        r.inferred["formulation.integration"] = \
            "implied by the element formulation (not user input)"

    # --- solution field
    if not r.fields.solution:
        r.missing.append("fields.solution: the converged solution field U "
                         "(e.g. U.npy, exported from the ODB)")

    # --- material evaluator
    if not (r.material.file or r.material.type or r.material.properties):
        needs_mat = any(b.startswith("solid_") for b in r.formulation.backends.values())
        if needs_mat:
            r.missing.append(
                "material: a material evaluator (UMAT source, or type + properties) "
                "for the solid elements")

    # --- constraints / loads come from the mesh when it is an .inp
    if r.mesh.format == "abaqus_inp":
        nb = len(getattr(model, "boundaries", []) or [])
        nl = len(getattr(model, "cloads", []) or [])
        if nb:
            r.inferred["constraints"] = "%d *Boundary block(s) read from the mesh" % nb
        if nl:
            r.inferred["stimuli.loads"] = "%d *Cload(s) read from the mesh" % nl

    # --- sensitivity defaults
    r.sensitivity.setdefault("order", 1)
    r.sensitivity.setdefault("backend", "otilib")
    r.inferred.setdefault("sensitivity.order", "default 1")
    r.inferred.setdefault("sensitivity.backend", "default otilib")
    r.inferred.setdefault("output.dir", "default '%s'" % r.output_dir)

    if not r.parameters:
        r.missing.append("parameters: which parameters to differentiate "
                         "(e.g. 'MAT.E', 'spring.k')")
    return r


def _load_model(r: Recipe):
    """Parse the mesh into a neutral core Model."""
    path = r.path(r.mesh.file)
    if r.mesh.format == "neutral_json":
        from residual_core.io.neutral_model_io import load_model
        return load_model(path)
    if r.mesh.format == "abaqus_inp":
        from residual_core.io import abaqus_inp_parser as P
        from residual_core.core.model import from_abaqus
        return from_abaqus(P.parse_inp(path))
    raise ValueError("unknown mesh format %r (use .inp or .json)" % r.mesh.format)


# --------------------------------------------------------------------------- #
# Capability gating — the honest part
# --------------------------------------------------------------------------- #
def sensitivity_capability(r: Recipe) -> Dict[str, Any]:
    """Can we actually OTI-differentiate this assembled residual?

    Assembling R and DIFFERENTIATING R are different capabilities. Most backends
    can do the first and not the second: their kernels use numpy float arrays,
    which cannot carry an OTI number. We refuse to pretend otherwise.
    """
    from residual_core.formulations import default_formulations
    reg = default_formulations()
    ok, bad = [], []
    for et, key in sorted(r.formulation.backends.items()):
        f = reg.get(key)
        (ok if getattr(f, "oti_differentiable", False) else bad).append((et, key))

    out = {"can_assemble": bool(r.formulation.backends) and not r.blockers,
           "oti_differentiable": bool(ok) and not bad,
           "oti_ready": [("%s -> %s" % (e, k)) for e, k in ok],
           "oti_blocked": [("%s -> %s" % (e, k)) for e, k in bad]}
    if bad:
        out["reason"] = (
            "these element formulations cannot be evaluated with OTI scalars: %s. "
            "Their kernels use numpy FLOAT arrays (e.g. core/voigt.py::isotropic_D "
            "raises on an OTI number), so parameters cannot be overloaded through "
            "them. You can still ASSEMBLE and verify R. For sensitivities today, "
            "use the black-box path (your solver returns the residual "
            "coefficients), or an OTI-transformed UMAT. Making the solid kernels "
            "OTI-safe is a known engineering task, not a physics limit."
            % ", ".join("%s (%s)" % (e, k) for e, k in bad))
    return out


def describe(r: Recipe) -> str:
    """Human-readable: what you gave, what we inferred, what is missing."""
    L = ["Residual assembly recipe: %s" % r.name, "-" * 46]
    L.append("Mesh:")
    L.append("  file            : %s" % (r.mesh.file or "-- MISSING --"))
    if r.mesh.n_nodes is not None:
        L.append("  nodes/elements  : %d / %d" % (r.mesh.n_nodes, r.mesh.n_elements))
        L.append("  element types   : %s" % ", ".join(r.mesh.element_types))
    L.append("Ingredients:")
    L.append("  solution field  : %s" % (r.fields.solution or "-- MISSING --"))
    L.append("  material        : %s" % (r.material.file or r.material.type or
                                         (", ".join(r.material.properties) or "-- MISSING --")))
    L.append("  element fields  : %s" % (r.fields.element_fields or "(not given)"))
    L.append("  parameters      : %s" % (", ".join(r.parameters) or "-- MISSING --"))

    if r.inferred:
        L.append("")
        L.append("Inferred for you (you did not have to type these):")
        for k, v in sorted(r.inferred.items()):
            L.append("  %-22s %s" % (k, v))

    cap = sensitivity_capability(r)
    L.append("")
    L.append("Capability:")
    L.append("  assemble R           : %s" % ("yes" if cap["can_assemble"] else "no"))
    L.append("  OTI-differentiate R  : %s" % ("yes" if cap["oti_differentiable"] else "NO"))
    if cap.get("oti_blocked"):
        L.append("     blocked: %s" % ", ".join(cap["oti_blocked"]))

    if r.blockers:
        L.append("")
        L.append("Blockers:")
        for b in r.blockers:
            L.append("  [blocked] %s" % b)
    if r.missing:
        L.append("")
        L.append("Minimum missing input:")
        L.append("  [missing] %s" % r.missing[0])
        for m in r.missing[1:]:
            L.append("  [missing] %s" % m)
    if not r.missing and not r.blockers:
        L.append("")
        L.append("Ready to assemble.")
    return "\n".join(L)
