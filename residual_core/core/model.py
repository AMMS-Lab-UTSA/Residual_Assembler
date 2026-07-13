"""Formulation-agnostic model container.

`Model` holds only mesh + bindings + boundary data. It knows which *formulation*
and which *material* each element is bound to (by registry key), but nothing
about what those formulations/materials actually do. The assembler consumes this
without any crystal-plasticity (or any other physics) assumptions.

`from_abaqus` converts the Abaqus-specific parse result (io/abaqus_inp_parser's
AbaqusModel) into this neutral form, choosing a formulation per element via a
caller-supplied policy (default: C3D8 -> solid_c3d8_finite_strain).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np


@dataclass
class Element:
    eid: int
    etype: str
    connectivity: List[int]


@dataclass
class Model:
    nodes: Dict[int, Tuple[float, float, float]] = field(default_factory=dict)
    elements: Dict[int, Element] = field(default_factory=dict)
    node_sets: Dict[str, List[int]] = field(default_factory=dict)
    element_sets: Dict[str, List[int]] = field(default_factory=dict)

    # bindings (registry keys, not objects -> the assembler resolves them)
    element_formulation: Dict[int, str] = field(default_factory=dict)
    element_material: Dict[int, str] = field(default_factory=dict)

    # named materials -> a MaterialBinding-like object or a raw dict of props
    materials: Dict[str, Any] = field(default_factory=dict)

    # boundary / load / constraint data (generic)
    boundaries: List[Any] = field(default_factory=list)
    cloads: List[Any] = field(default_factory=list)
    equations: List[Any] = field(default_factory=list)

    ndim: int = 3
    meta: Dict[str, Any] = field(default_factory=dict)

    def coords_of(self, connectivity) -> np.ndarray:
        return np.array([self.nodes[nid] for nid in connectivity], dtype=float)

    def elements_of_type(self, etype):
        return [e for e in self.elements.values() if e.etype.upper() == etype.upper()]


def default_formulation_policy(etype: str) -> Optional[str]:
    """Map an element type to a formulation registry key. Returns None for
    element types with no registered continuum backend (kept, but not assembled
    by the default solid path)."""
    et = etype.upper()
    if et == "C3D8":
        return "solid_c3d8_finite_strain"
    # C3D8R / C3D20R / C3D4 / U1 / cohesive: no default continuum backend yet
    return None


def from_abaqus(am, formulation_policy: Callable[[str], Optional[str]] = None,
                material_bindings: Optional[Dict[str, Any]] = None) -> Model:
    """Convert an AbaqusModel (io/abaqus_inp_parser) into a neutral Model.

    formulation_policy : etype -> formulation key (default: C3D8 finite strain).
    material_bindings  : optional {material_name -> MaterialBinding}; if omitted,
                         the raw parsed Material objects are stored and a binding
                         must be attached before a Mode-2 assembly.
    """
    policy = formulation_policy or default_formulation_policy
    m = Model()
    m.nodes = dict(am.nodes)
    for eid, el in am.elements.items():
        m.elements[eid] = Element(eid, el.etype, list(el.connectivity))
        fk = policy(el.etype)
        if fk is not None:
            m.element_formulation[eid] = fk
    m.node_sets = {k: list(v) for k, v in am.node_sets.items()}
    m.element_sets = {k: list(v) for k, v in am.element_sets.items()}
    m.element_material = dict(am.element_material)
    m.materials = dict(material_bindings) if material_bindings else dict(am.materials)
    m.boundaries = list(am.boundaries)
    m.cloads = list(am.cloads)
    m.equations = list(am.equations)
    m.meta = {"source": "abaqus_inp", "n_c3d8": len(m.elements_of_type("C3D8"))}
    return m
