"""Model inspector / diagnostics — the auto-detection UX.

Given a parsed model, decide (using the registries only — never physics) which
element types and materials can be assembled, by which backend, in which modes,
and report the *minimum* missing user inputs. This is the engine behind
``resasm inspect`` / ``resasm doctor``.

Example rendered output::

    Model inspection summary
    ------------------------
    Elements:
      - B31 beam elements: supported by beam2 backend
      - S4 shell elements: backend missing
      - C3D8 solid elements: supported by solid_c3d8_finite_strain backend

    Materials:
      - ElasticMaterial1: built-in elastic adapter available
      - UserMaterial2: UMAT adapter required; source file not provided

    Required user inputs:
      1. Provide UMAT source or stress field export for UserMaterial2.
      2. Shell backend for S4 is not implemented yet.

    Possible modes:
      - stress-driven residual: available if exported fields are provided (C3D8)
      - material replay: available for C3D8
      - UEL-direct: not applicable
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class ElementReport:
    etype: str
    count: int
    backends: List[str] = field(default_factory=list)   # names that can handle it
    implemented: bool = True     # a matched backend that is actually runnable
    status: str = "unspecified"                         # verification_status
    limitations: List[str] = field(default_factory=list)
    modes: List[str] = field(default_factory=list)      # available assembly modes
    next_step: str = ""                                 # minimum next step if blocked

    @property
    def supported(self) -> bool:
        return bool(self.backends) and self.implemented

    @property
    def selected_backend(self) -> Optional[str]:
        return self.backends[0] if self.backends else None


@dataclass
class MaterialReport:
    name: str
    is_user_material: bool = False
    backend: Optional[str] = None       # resolved material backend name
    source_provided: bool = False       # UMAT source / adapter attached?

    @property
    def status(self) -> str:
        if not self.is_user_material:
            return "built-in material adapter available" if self.backend \
                else "built-in material (no adapter matched)"
        if self.source_provided:
            return "UMAT adapter available; source attached"
        return "UMAT adapter required; source file not provided"


@dataclass
class InspectionReport:
    elements: List[ElementReport] = field(default_factory=list)
    materials: List[MaterialReport] = field(default_factory=list)
    required_inputs: List[str] = field(default_factory=list)
    modes: Dict[str, str] = field(default_factory=dict)   # mode -> availability text

    def render(self) -> str:
        L = ["Model inspection summary", "------------------------", "Elements:"]
        for e in self.elements:
            if e.backends and e.implemented:
                L.append("  - %s elements (x%d): supported by %s backend"
                         % (e.etype, e.count, e.backends[0]))
            elif e.backends and not e.implemented:
                L.append("  - %s elements (x%d): recognised (%s) but backend not "
                         "implemented yet" % (e.etype, e.count, e.backends[0]))
            else:
                L.append("  - %s elements (x%d): backend missing" % (e.etype, e.count))
        L.append("")
        L.append("Materials:")
        if not self.materials:
            L.append("  (none)")
        for m in self.materials:
            L.append("  - %s: %s" % (m.name, m.status))
        L.append("")
        L.append("Required user inputs:")
        if not self.required_inputs:
            L.append("  (none — model is ready to assemble)")
        for i, req in enumerate(self.required_inputs, 1):
            L.append("  %d. %s" % (i, req))
        L.append("")
        L.append("Possible modes:")
        for mode, txt in self.modes.items():
            L.append("  - %s: %s" % (mode, txt))
        return "\n".join(L)

    def render_backends(self) -> str:
        """Per-element-type backend selection detail (registry audit UX)."""
        L = ["Backend selection by element type",
             "---------------------------------"]
        for e in self.elements:
            L.append("Element type %s (x%d):" % (e.etype, e.count))
            L.append("  selected backend: %s" % (e.selected_backend or "(none)"))
            L.append("  status: %s" % e.status)
            if e.limitations:
                L.append("  limitations: %s" % "; ".join(e.limitations))
            L.append("  available modes: %s"
                     % (", ".join(e.modes) if e.modes else "none"))
            if e.next_step:
                L.append("  minimum next step: %s" % e.next_step)
        return "\n".join(L)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "elements": [{"etype": e.etype, "count": e.count,
                          "backends": e.backends, "supported": e.supported,
                          "status": e.status, "limitations": e.limitations,
                          "modes": e.modes, "next_step": e.next_step}
                         for e in self.elements],
            "materials": [{"name": m.name, "user_material": m.is_user_material,
                           "backend": m.backend, "status": m.status}
                          for m in self.materials],
            "required_inputs": list(self.required_inputs),
            "modes": dict(self.modes),
        }


def _iter_elements(model):
    """Yield (etype, ) for each element of a neutral Model or an AbaqusModel."""
    for el in model.elements.values():
        yield el.etype


def _material_items(model):
    """Yield (name, obj) for each material; obj may be a parser Material, a
    MaterialBinding, or a raw dict — duck-typed."""
    mats = getattr(model, "materials", {}) or {}
    for name, obj in mats.items():
        yield name, obj


def inspect_model(model, formulation_registry, material_registry=None,
                  attached_subroutine: bool = False) -> InspectionReport:
    """Build an InspectionReport from a model + the backend registries.

    ``formulation_registry`` / ``material_registry`` are ``core.registry.Registry``
    instances (or any object exposing ``find_for_element`` / ``get`` / ``specs``).
    ``attached_subroutine`` marks that a UMAT/UEL source has been provided.
    """
    rep = InspectionReport()

    # ---- elements --------------------------------------------------------
    counts = Counter(_iter_elements(model))
    for etype, n in sorted(counts.items()):
        backends = formulation_registry.find_for_element(etype)
        implemented = True
        status = "no backend"
        limitations: List[str] = []
        modes: List[str] = []
        next_step = ""
        if backends:
            spec = getattr(formulation_registry.get(backends[0]), "spec", None)
            implemented = bool(spec and spec.supported_modes)
            if spec:
                status = spec.verification_status or "unspecified"
                limitations = list(spec.limitations)
            # available modes = union across ALL backends matching this etype
            seen = set()
            for b in backends:
                bspec = getattr(formulation_registry.get(b), "spec", None)
                for md in (bspec.supported_modes if bspec else ()):
                    if md not in seen:
                        seen.add(md)
                        modes.append(md)
            if not implemented:
                next_step = ("implement/register a backend for %s, or use "
                             "UEL-direct mode with a callable element routine"
                             % etype)
        else:
            next_step = ("register a formulation backend that declares element "
                         "type %s" % etype)
        rep.elements.append(ElementReport(
            etype, n, backends, implemented, status=status,
            limitations=limitations, modes=modes, next_step=next_step))

    # ---- materials -------------------------------------------------------
    for name, obj in _material_items(model):
        is_user = bool(getattr(obj, "user_material", False))
        backend = None
        if material_registry is not None:
            if is_user and "umat" in material_registry:
                backend = "umat"
            elif not is_user and "isotropic_elastic" in material_registry:
                backend = "isotropic_elastic"
        rep.materials.append(MaterialReport(
            name=name, is_user_material=is_user, backend=backend,
            source_provided=attached_subroutine))

    # ---- required user inputs (minimum, not a generic checklist) ---------
    for e in rep.elements:
        if not e.backends:
            rep.required_inputs.append(
                "No backend for element type %s; register a formulation that "
                "declares it." % e.etype)
        elif not e.implemented:
            rep.required_inputs.append(
                "Backend for %s (%s) is a contract placeholder — implementation "
                "not complete." % (e.etype, e.backends[0]))
    for m in rep.materials:
        if m.is_user_material and not m.source_provided:
            rep.required_inputs.append(
                "Provide UMAT source or an exported field for material '%s' "
                "(or run stress-driven mode)." % m.name)

    # ---- possible modes --------------------------------------------------
    def _etypes_supporting(mode):
        out = []
        for e in rep.elements:
            for b in e.backends:
                spec = getattr(formulation_registry.get(b), "spec", None)
                if spec and mode in spec.supported_modes:
                    out.append(e.etype)
                    break
        return sorted(set(out))

    sd = _etypes_supporting("stress-driven")
    mr = _etypes_supporting("material-replay")
    fm = _etypes_supporting("formulation")
    uel = [e.etype for e in rep.elements
           if any(formulation_registry.get(b) is not None
                  and getattr(formulation_registry.get(b), "name", "") == "uel_direct"
                  for b in e.backends)]

    rep.modes["formulation (self-contained)"] = (
        "available for %s (needs section/constitutive properties)"
        % ", ".join(fm)) if fm else "not applicable"
    rep.modes["stress-driven residual"] = (
        "available if exported stress/resultant fields are provided (%s)"
        % ", ".join(sd)) if sd else "not applicable"
    rep.modes["material replay"] = (
        "available for %s" % ", ".join(mr)) if mr else "not applicable"
    rep.modes["UEL-direct"] = (
        "available for %s" % ", ".join(sorted(set(uel)))) if uel else "not applicable"
    return rep
