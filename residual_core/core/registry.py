"""Generic registry + backend capability descriptor.

The core is model-agnostic: it does not import concrete formulations or
materials. Instead each backend *declares* what it can do (a ``BackendSpec``)
and registers itself. The model inspector, requirements engine and CLI then reason
purely over these declarations — never over physics.

This module is intentionally free of any FEM, mechanics, or material knowledge.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Tuple


@dataclass
class BackendSpec:
    """Machine-readable declaration of what a backend supports.

    Every formulation and material advertises one of these so the inspector can
    decide — without knowing any physics — which elements/materials can be
    assembled, in which modes, and what inputs are still missing.
    """
    backend_name: str
    kind: str                                   # "formulation" | "material"
    supported_element_types: Tuple[str, ...] = ()
    dof_types: Tuple[str, ...] = ()
    supported_modes: Tuple[str, ...] = ()
    required_inputs: Tuple[str, ...] = ()
    optional_inputs: Tuple[str, ...] = ()
    verification_tests: Tuple[str, ...] = ()
    limitations: Tuple[str, ...] = ()
    notes: str = ""

    # ---- richer, per-mode capability declaration (audit item 3) ----------
    # Inputs can differ per mode; declare them explicitly so the requirements
    # engine and the CLI never have to guess. Falls back to the flat
    # required_inputs / optional_inputs when a mode is absent here.
    required_inputs_by_mode: Dict[str, Tuple[str, ...]] = field(default_factory=dict)
    optional_inputs_by_mode: Dict[str, Tuple[str, ...]] = field(default_factory=dict)
    material_interface_needed: bool = False     # does it call a Material.evaluate?
    state_requirements: str = "none"            # none | history-dependent | ...
    tangent_support: str = "none"               # analytic | fd | none
    verification_status: str = "unspecified"    # verified | simple | contract-only | skeleton

    def inputs_for_mode(self, mode: str) -> Tuple[str, ...]:
        """Required inputs for a specific mode (falls back to the flat list)."""
        return tuple(self.required_inputs_by_mode.get(mode, self.required_inputs))

    def optional_for_mode(self, mode: str) -> Tuple[str, ...]:
        return tuple(self.optional_inputs_by_mode.get(mode, self.optional_inputs))

    def matches_element(self, etype: str) -> bool:
        """True if this backend handles ``etype``. Supports ``*`` suffix
        wildcards (e.g. ``U*`` matches any Abaqus user element)."""
        et = etype.upper()
        for pat in self.supported_element_types:
            p = pat.upper()
            if p.endswith("*") and et.startswith(p[:-1]):
                return True
            if p == et:
                return True
        return False

    def as_dict(self) -> Dict[str, Any]:
        return {
            "backend_name": self.backend_name,
            "kind": self.kind,
            "supported_element_types": list(self.supported_element_types),
            "dof_types": list(self.dof_types),
            "supported_modes": list(self.supported_modes),
            "required_inputs": list(self.required_inputs),
            "optional_inputs": list(self.optional_inputs),
            "verification_tests": list(self.verification_tests),
            "limitations": list(self.limitations),
            "notes": self.notes,
            "required_inputs_by_mode": {k: list(v) for k, v in
                                        self.required_inputs_by_mode.items()},
            "optional_inputs_by_mode": {k: list(v) for k, v in
                                        self.optional_inputs_by_mode.items()},
            "material_interface_needed": self.material_interface_needed,
            "state_requirements": self.state_requirements,
            "tangent_support": self.tangent_support,
            "verification_status": self.verification_status,
        }


class Registry:
    """A name -> instance store whose entries expose a ``.spec`` BackendSpec.

    Kept deliberately tiny and physics-free. Both the formulation registry and
    the material registry are instances of this class.
    """

    def __init__(self, kind: str):
        self.kind = kind
        self._items: Dict[str, Any] = {}

    def register(self, instance: Any, name: Optional[str] = None) -> Any:
        key = name or getattr(instance, "name", None) or getattr(
            getattr(instance, "spec", None), "backend_name", None)
        if not key:
            raise ValueError("cannot register a %s backend without a name" % self.kind)
        self._items[key] = instance
        return instance

    def get(self, name: str) -> Any:
        return self._items.get(name)

    def __contains__(self, name: str) -> bool:
        return name in self._items

    def __getitem__(self, name: str) -> Any:
        return self._items[name]

    def names(self) -> List[str]:
        return sorted(self._items)

    def instances(self) -> Iterable[Any]:
        return list(self._items.values())

    def specs(self) -> List[BackendSpec]:
        out = []
        for inst in self._items.values():
            spec = getattr(inst, "spec", None)
            if spec is not None:
                out.append(spec)
        return out

    def as_dict(self) -> Dict[str, Any]:
        return {name: (getattr(inst, "spec", None).as_dict()
                       if getattr(inst, "spec", None) else {})
                for name, inst in self._items.items()}

    def find_for_element(self, etype: str) -> List[str]:
        """Names of registered backends that can handle ``etype``."""
        out = []
        for name, inst in self._items.items():
            spec = getattr(inst, "spec", None)
            if spec is not None and spec.matches_element(etype):
                out.append(name)
        return sorted(out)
