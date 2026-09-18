"""Inputs of the history replay: the Abaqus deck subset and recorded fields.

``read_history_model`` reads an ``Analysis.inp`` into a :class:`HistoryModel`
and refuses anything the replay cannot reproduce: every keyword outside the
supported subset, amplitudes, ``OP=NEW``, NLGEOM, element types other than
C3D8, more than one material or step. It never guesses a missing feature.

Supported (one static step, small strain):

- ``*Node``, ``*Element, type=C3D8``, ``*Nset``/``*Elset`` (incl. GENERATE),
  one ``*Part``/``*Instance`` without transformation, ``*Material`` with
  ``*User Material`` and ``*Depvar``, one ``*Solid Section`` covering every
  element;
- ``*Boundary`` in the model data: zero-valued (ENCASTRE/PINNED/XSYMM/...,
  or DOF ranges with value 0);
- ``*Boundary`` inside the step: any value, ramped linearly over the step
  from the value in force at the step start (Abaqus's default ramp amplitude
  for a static step), ``OP=MOD``;
- ``*Cload`` inside the step: ramped from zero, ``OP=MOD``;
- ``*Static`` (the increments actually taken come from the ODB record),
  ``*Controls`` (solver controls only), output requests.

``load_recorded_fields`` reads the per-frame export of the ODB (``.npz``,
written by :mod:`residual_core.replay.odb_export_npz` or the equivalent
cantilever exporter) and checks it against the deck: node and element labels,
coordinates (to single precision), connectivity, a virgin frame 0, strictly
increasing times that end at the step period, and prescribed displacements
equal to the deck's ramped values.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

import numpy as np

from ..io.abaqus_inp_parser import _iter_blocks, parse_inp

SUPPORTED_KEYWORDS = {
    "heading": set(), "preprint": {"echo", "model", "history", "contact"},
    "node": {"nset"}, "element": {"type", "elset"},
    "nset": {"nset", "generate", "instance", "internal", "unsorted"},
    "elset": {"elset", "generate", "instance", "internal", "unsorted"},
    "part": {"name"}, "end part": set(), "assembly": {"name"}, "end assembly": set(),
    "instance": {"name", "part"}, "end instance": set(),
    "material": {"name"}, "user material": {"constants", "type"}, "depvar": set(),
    "solid section": {"elset", "material"}, "boundary": {"op"}, "cload": {"op"},
    "step": {"name", "nlgeom", "inc"}, "static": {"direct"}, "controls": {"parameters", "reset"},
    "end step": set(), "output": {"field", "history", "frequency", "variable"},
    "node output": {"nset", "frequency"}, "element output": {"elset", "position", "directions", "frequency"},
}
PRECISION_FLOOR = float(np.finfo(np.float32).eps)


class UnsupportedFeature(ValueError):
    """The deck or record uses a feature the history replay does not reproduce."""


@dataclass
class HistoryModel:
    node_ids: np.ndarray
    coords: np.ndarray
    element_ids: np.ndarray
    connectivity: np.ndarray                  # (ne, 8) zero-based node positions
    node_sets: Dict[str, List[int]]
    element_sets: Dict[str, List[int]]
    props: np.ndarray
    depvar: int
    material_name: str
    step_name: str
    time_period: float
    prescribed_dofs: np.ndarray               # (nc,) global DOF indices, sorted
    prescribed_start: np.ndarray              # value at the step start
    prescribed_end: np.ndarray                # value at the step end
    load_end: np.ndarray                      # (ndof,) concentrated loads at the step end
    source: str = ""
    notes: List[str] = field(default_factory=list)

    @property
    def ndof(self) -> int:
        return 3 * len(self.node_ids)

    def prescribed_at(self, time: float) -> np.ndarray:
        """Prescribed DOF values at step time ``time`` (Abaqus static ramp)."""
        factor = time / self.time_period
        return self.prescribed_start + (self.prescribed_end - self.prescribed_start) * factor

    def load_at(self, time: float) -> np.ndarray:
        return self.load_end * (time / self.time_period)

    def node_set(self, name: str) -> List[int]:
        key = name.upper()
        if key not in self.node_sets:
            raise ValueError("unknown node set %r (available: %s)" % (name, sorted(self.node_sets)))
        return self.node_sets[key]

    def element_set(self, name: str) -> List[int]:
        key = name.upper()
        if key not in self.element_sets:
            raise ValueError("unknown element set %r (available: %s)" % (name, sorted(self.element_sets)))
        return self.element_sets[key]


def _check_keywords(path: Path) -> None:
    counts: Dict[str, int] = {}
    for keyword, spelling, options, data in _iter_blocks(path.read_text().splitlines()):
        counts[keyword] = counts.get(keyword, 0) + 1
        if keyword not in SUPPORTED_KEYWORDS:
            raise UnsupportedFeature("unsupported INP keyword *%s" % spelling)
        unknown = set(options) - SUPPORTED_KEYWORDS[keyword]
        if unknown:
            raise UnsupportedFeature("unsupported *%s options: %s" % (spelling, sorted(unknown)))
        if keyword == "element" and str(options.get("type", "")).upper() != "C3D8":
            raise UnsupportedFeature("element type %s: the history replay supports C3D8 "
                                     "(selective-reduced B-bar) only" % options.get("type"))
        if keyword == "user material" and str(options.get("type", "MECHANICAL")).upper() != "MECHANICAL":
            raise UnsupportedFeature("only mechanical user materials are supported")
        if keyword in ("boundary", "cload") and str(options.get("op", "MOD")).upper() != "MOD":
            raise UnsupportedFeature("*%s OP=NEW is not supported" % spelling)
        if keyword == "solid section" and any(line.replace(",", "").strip() for line in data):
            raise UnsupportedFeature("*Solid Section data lines (thickness etc.) are not supported")
    for keyword in ("step", "static", "end step", "material", "user material", "solid section"):
        if counts.get(keyword, 0) != 1:
            raise UnsupportedFeature("exactly one *%s block is required (found %d)"
                                     % (keyword, counts.get(keyword, 0)))
    if counts.get("part", 0) > 1 or counts.get("instance", 0) > 1:
        raise UnsupportedFeature("only one part/instance is supported")


def _targets(model, target) -> List[int]:
    if isinstance(target, int):
        if target not in model.nodes:
            raise ValueError("boundary/load node %d does not exist" % target)
        return [target]
    for name, members in model.node_sets.items():
        if name.upper() == str(target).upper():
            if not members:
                raise ValueError("node set %s is empty" % target)
            return list(members)
    raise ValueError("unknown node set %r" % target)


def read_history_model(path) -> HistoryModel:
    path = Path(path)
    if not path.is_file():
        raise ValueError("model deck not found: %s" % path)
    _check_keywords(path)
    model = parse_inp(str(path))
    if model.warnings:
        raise UnsupportedFeature("INP parser diagnostics: " + "; ".join(model.warnings))
    # parse_inp lists every keyword it does not model (output requests,
    # *Controls, *Heading ...); _check_keywords above already refused anything
    # outside SUPPORTED_KEYWORDS, so that list carries no further information.
    step = model.steps[0]
    if step.nlgeom:
        raise UnsupportedFeature("NLGEOM=YES: the history replay is small strain only")
    if step.procedure != "static":
        raise UnsupportedFeature("the step must be *Static")
    if not step.time_period or not np.isfinite(step.time_period) or step.time_period <= 0:
        raise ValueError("*Static needs a positive step time period")
    if len(model.materials) != 1:
        raise UnsupportedFeature("exactly one material is supported")
    material = next(iter(model.materials.values()))
    if not material.user_material or not material.constants:
        raise UnsupportedFeature("the material must be a *User Material with constants")
    if set(model.element_material) != set(model.elements):
        raise UnsupportedFeature("every element must be assigned the user material by one solid section")
    node_ids = np.array(sorted(model.nodes), dtype=int)
    position = {int(node): index for index, node in enumerate(node_ids)}
    coords = np.array([model.nodes[node] for node in node_ids], dtype=float)
    element_ids = np.array(sorted(model.elements), dtype=int)
    connectivity = np.array([[position[node] for node in model.elements[eid].connectivity]
                             for eid in element_ids], dtype=int)
    if connectivity.shape[1:] != (8,):
        raise UnsupportedFeature("C3D8 elements need 8 nodes")

    values: Dict[int, List[float]] = {}          # dof -> [start, end]
    for boundary in model.boundaries:
        if boundary.amplitude:
            raise UnsupportedFeature("boundary amplitudes are not supported (ramp only)")
        if boundary.step not in (0, 1):
            raise UnsupportedFeature("boundary in a second step")
        first, last = boundary.dof_start, min(boundary.dof_end, 3)
        if first < 1 or first > 3 or last < first:
            raise UnsupportedFeature("boundary DOFs %d-%d: solid elements have DOFs 1-3"
                                     % (boundary.dof_start, boundary.dof_end))
        if boundary.step == 0 and boundary.value != 0:
            raise UnsupportedFeature("nonzero model-data boundary values are not supported; "
                                     "prescribe them inside the step")
        for node in _targets(model, boundary.target):
            for dof in range(first, last + 1):
                index = 3 * position[node] + dof - 1
                start_end = values.setdefault(index, [0.0, 0.0])
                if boundary.step == 1:
                    start_end[1] = float(boundary.value)
    for index, (start, end) in values.items():
        if not np.isfinite(end):
            raise ValueError("nonfinite boundary value at DOF %d" % index)
    prescribed = np.array(sorted(values), dtype=int)
    load = np.zeros(3 * len(node_ids))
    for cload in model.cloads:
        if cload.amplitude:
            raise UnsupportedFeature("load amplitudes are not supported (ramp only)")
        if cload.step != 1:
            raise UnsupportedFeature("concentrated loads must be defined inside the step")
        if cload.dof not in (1, 2, 3) or not np.isfinite(cload.value):
            raise UnsupportedFeature("concentrated loads need DOF 1-3 and a finite value")
        for node in _targets(model, cload.target):
            load[3 * position[node] + cload.dof - 1] += cload.value
    if np.any(load[prescribed] != 0):
        raise UnsupportedFeature("a concentrated load acts on a prescribed DOF")
    if not len(prescribed):
        raise ValueError("no boundary conditions: the stiffness is singular")
    return HistoryModel(
        node_ids=node_ids, coords=coords, element_ids=element_ids, connectivity=connectivity,
        node_sets={name.upper(): sorted(set(int(n) for n in members)) for name, members in model.node_sets.items()},
        element_sets={name.upper(): sorted(set(int(e) for e in members)) for name, members in model.element_sets.items()},
        props=np.array(material.constants, dtype=float), depvar=int(material.depvar or 0),
        material_name=material.name, step_name=step.name, time_period=float(step.time_period),
        prescribed_dofs=prescribed,
        prescribed_start=np.array([values[i][0] for i in prescribed]),
        prescribed_end=np.array([values[i][1] for i in prescribed]),
        load_end=load, source=str(path))


@dataclass
class RecordedFields:
    time: np.ndarray            # (nf,) step time, frame 0 = virgin
    U: np.ndarray               # (nf, nn, 3)
    RF: np.ndarray              # (nf, nn, 3)
    S: np.ndarray               # (nf, ne, 8, 6)
    SDV: np.ndarray             # (nf, ne, 8, nsdv)
    source: str
    precision: str = "float32"  # ODB field output precision
    CF: np.ndarray = None

    @property
    def increments(self) -> int:
        return len(self.time) - 1


def load_recorded_fields(path, model: HistoryModel) -> RecordedFields:
    """Load an ODB field export and check it matches the deck exactly."""
    path = Path(path)
    data = np.load(path, allow_pickle=False)
    required = {"time", "node_labels", "coords", "elem_labels", "conn", "U", "RF", "S", "SDV"}
    missing = required - set(data.files)
    if missing:
        raise ValueError("field export %s lacks %s" % (path, sorted(missing)))
    node_labels = np.asarray(data["node_labels"], dtype=int)
    elem_labels = np.asarray(data["elem_labels"], dtype=int)
    if not np.array_equal(np.sort(node_labels), model.node_ids):
        raise ValueError("INP/ODB node labels differ")
    if not np.array_equal(np.sort(elem_labels), model.element_ids):
        raise ValueError("INP/ODB element labels differ")
    node_order = np.argsort(node_labels)
    elem_order = np.argsort(elem_labels)
    coords = np.asarray(data["coords"], dtype=float)[node_order]
    scale = max(np.abs(model.coords).max(), 1.0)
    if np.abs(coords - model.coords).max() > 4 * PRECISION_FLOOR * scale:
        raise ValueError("INP/ODB node coordinates differ beyond single precision")
    conn = np.asarray(data["conn"], dtype=int)[elem_order]
    if not np.array_equal(conn, model.node_ids[model.connectivity]):
        raise ValueError("INP/ODB connectivity differs")
    time = np.asarray(data["time"], dtype=float)
    if time.ndim != 1 or len(time) < 2 or time[0] != 0 or np.any(np.diff(time) <= 0):
        raise ValueError("the export needs a virgin frame 0 and strictly increasing frame times")
    if abs(time[-1] - model.time_period) > 1e-6 * model.time_period:
        raise ValueError("the last frame (t=%g) does not reach the step period %g"
                         % (time[-1], model.time_period))
    U = np.asarray(data["U"], dtype=float)[:, node_order]
    RF = np.asarray(data["RF"], dtype=float)[:, node_order]
    S = np.asarray(data["S"], dtype=float)[:, elem_order]
    SDV = np.asarray(data["SDV"], dtype=float)[:, elem_order]
    nf, nn, ne = len(time), len(model.node_ids), len(model.element_ids)
    for name, array, shape in (("U", U, (nf, nn, 3)), ("RF", RF, (nf, nn, 3)),
                               ("S", S, (nf, ne, 8, 6))):
        if array.shape != shape or not np.all(np.isfinite(array)):
            raise ValueError("export %s must be finite with shape %s, got %s" % (name, shape, array.shape))
    if SDV.shape[:3] != (nf, ne, 8) or SDV.shape[3] != model.depvar:
        raise ValueError("export SDV has %d components; the deck declares *Depvar %d"
                         % (SDV.shape[3] if SDV.ndim == 4 else -1, model.depvar))
    if np.any(U[0]) or np.any(S[0]) or np.any(SDV[0]):
        raise UnsupportedFeature("nonzero initial state (frame 0) is not supported")
    CF = None
    if "CF" in data.files:
        CF = np.asarray(data["CF"], dtype=float)[:, node_order]
        expected = np.array([model.load_at(t) for t in time]).reshape(nf, nn, 3)
        if np.abs(CF - expected).max() > 4 * PRECISION_FLOOR * max(np.abs(expected).max(), 1e-30) + 1e-30:
            raise ValueError("ODB concentrated loads (CF) differ from the deck's ramped *Cload")
    elif np.any(model.load_end):
        raise ValueError("the deck has concentrated loads; the export must include CF to check them")
    flat = U.reshape(nf, -1)
    for frame in range(nf):
        target = model.prescribed_at(time[frame])
        recorded = flat[frame, model.prescribed_dofs]
        bound = 4 * PRECISION_FLOOR * max(np.abs(target).max(), 1e-30) + 1e-30
        if np.abs(recorded - target).max() > bound:
            raise ValueError("frame %d: recorded displacements differ from the deck's prescribed "
                             "values by %.3e (> %.3e)" % (frame, np.abs(recorded - target).max(), bound))
    return RecordedFields(time=time, U=U, RF=RF, S=S, SDV=SDV, source=str(path), CF=CF)
