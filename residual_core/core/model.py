"""Formulation-agnostic model container.

`Model` holds only mesh + bindings + boundary data. It knows which *formulation*
and which *material* each element is bound to (by registry key), but nothing
about what those formulations/materials actually do. The assembler consumes this
without any crystal-plasticity (or any other physics) assumptions.

`from_abaqus` converts the Abaqus-specific parse result (io/abaqus_inp_parser's
AbaqusModel) into this neutral form, choosing a formulation per element via a
caller-supplied policy (default: C3D8 -> solid_c3d8_finite_strain). A deck
keyword that can change equilibrium but that this model does not apply
(``*Dsload``, ``*Dload``, ``*Equation``, ``*Amplitude`` and ``AMPLITUDE=``
references, or any keyword the reader leaves unread) is never dropped
silently: ``from_abaqus`` names each one in ``Model.unapplied_keywords`` and in
a :class:`DeckKeywordNotApplied` warning.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np


class DeckKeywordNotApplied(UserWarning):
    """A deck keyword that can change equilibrium is present but not applied.

    Issued by :func:`from_abaqus` once per converted deck; the same notes are
    kept on ``Model.unapplied_keywords`` for ``resasm inspect``.
    """


#: Keywords that only request output, echo input, write restart data, tune
#: the solver or name a surface (whatever uses the surface is named itself).
#: They never enter R = F_int - F_ext, so leaving them unread changes nothing.
#: Every other keyword the reader does not apply is named.
RESIDUAL_NEUTRAL_KEYWORDS = frozenset({
    "heading", "preprint", "end step", "output", "node output", "element output",
    "energy output", "contact output", "node print", "el print", "energy print",
    "node file", "el file", "restart", "monitor", "controls", "file format",
    "surface"})


def unapplied_deck_keywords(am) -> List[str]:
    """What an Abaqus deck asks for that the general assembly does not apply.

    ``am`` is an ``io.abaqus_inp_parser.AbaqusModel``. One note per keyword,
    each starting with the keyword itself. Loads, constraints and amplitudes
    that change equilibrium, and every other keyword the reader left unread
    (apart from :data:`RESIDUAL_NEUTRAL_KEYWORDS`), are named here; the notes
    are empty for a deck whose every keyword is applied.
    """
    notes: List[str] = []
    dsloads = list(getattr(am, "dsloads", ()) or ())
    if dsloads:
        notes.append("*Dsload (%d line(s)): parsed but not applied; the residual "
                     "carries no distributed surface load" % len(dsloads))
    equations = list(getattr(am, "equations", ()) or ())
    if equations:
        notes.append("*Equation (%d constraint(s)): parsed but not applied; the "
                     "degrees of freedom it couples are treated as independent"
                     % len(equations))
    references: Dict[str, List[str]] = {}
    for kind, entries in (("*Boundary", getattr(am, "boundaries", ())),
                          ("*Cload", getattr(am, "cloads", ()))):
        for entry in entries or ():
            name = getattr(entry, "amplitude", None)
            if name:
                references.setdefault(str(name), []).append(kind)
    for name in sorted(references):
        kinds = references[name]
        notes.append("AMPLITUDE=%s (%s): not applied; these lines act at their full "
                     "magnitude" % (name, ", ".join(
                         "%d %s line(s)" % (kinds.count(kind), kind)
                         for kind in sorted(set(kinds)))))
    for keyword in getattr(am, "unsupported_keywords", ()) or ():
        normalized = " ".join(str(keyword).split()).lower()
        if normalized in RESIDUAL_NEUTRAL_KEYWORDS:
            continue
        if normalized == "dload":
            notes.append("*%s: present but not applied; the residual carries no "
                         "distributed or body load" % keyword)
        elif normalized == "amplitude":
            notes.append("*%s: present but not applied; loads and prescribed values "
                         "act at their full magnitude" % keyword)
        else:
            notes.append("*%s: present but not read, so not applied" % keyword)
    return notes


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
    #: The deck's ``*STEP`` blocks in order, when the source had any. Carried
    #: because a boundary condition without the step it was written in is a
    #: displacement with no time attached to it: a four-step verification deck
    #: read as one flat list applies the last step's displacements to the
    #: first step's increments.
    steps: List[Any] = field(default_factory=list)
    #: What the source deck asks for that this model does not apply (see
    #: :func:`unapplied_deck_keywords`); empty when every keyword is applied.
    unapplied_keywords: List[str] = field(default_factory=list)

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
    m.steps = list(getattr(am, "steps", ()) or ())
    m.unapplied_keywords = unapplied_deck_keywords(am)
    if m.unapplied_keywords:
        warnings.warn("the deck asks for what the residual assembly does not apply:\n  "
                      + "\n  ".join(m.unapplied_keywords),
                      DeckKeywordNotApplied, stacklevel=2)
    m.meta = {"source": "abaqus_inp", "n_c3d8": len(m.elements_of_type("C3D8"))}
    return m
