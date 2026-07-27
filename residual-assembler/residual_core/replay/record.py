"""Replay record: the sufficient primal analysis produced in Stage 1.

Defined by what the residual method mathematically NEEDS to reconstruct R, K and
R_,p at the converged solution -- not by what an ODB happens to contain. For a
small-strain elastic model the minimum is: mesh + connectivity + IP ordering,
the Dirichlet BCs and external loads, the material parameter values used, and
the converged nodal displacement (per required increment). Per-IP kinematics are
derived from the mesh + displacement; the production per-IP stress is carried
only to VALIDATE the replay.

For path-dependent models the record must additionally carry, per increment, the
exact material-point kinematic inputs and the initial state, so the replay can
march from the start propagating state derivatives. :func:`preflight_record`
reports exactly which required fields are missing for a given material package.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Optional

import numpy as np

from ..core.model import Element, Model

SCHEMA = "resasm_replay_record_v1"


class RecordError(Exception):
    """The replay record is unreadable or structurally invalid."""


@dataclass
class _BC:
    """Duck-typed boundary for core.constraints.partition."""
    target: int
    dof_start: int = 0
    dof_end: int = 0
    value: float = 0.0
    kind: str = "value"
    amplitude: Optional[str] = None


class ReplayRecord:
    def __init__(self, data: Mapping[str, Any]):
        self.raw = dict(data)

    @classmethod
    def load(cls, path: str) -> "ReplayRecord":
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
        except (OSError, ValueError) as exc:
            raise RecordError("cannot read replay record %r: %s" % (path, exc))
        return cls(data)

    # -- accessors -------------------------------------------------------- #
    @property
    def provenance(self) -> Dict[str, Any]:
        return dict(self.raw.get("provenance", {}))

    @property
    def kinematics(self) -> str:
        return str(self.raw.get("kinematics", ""))

    @property
    def props(self) -> List[float]:
        return [float(v) for v in self.raw.get("material", {}).get("props", [])]

    @property
    def increments(self) -> List[Dict[str, Any]]:
        return list(self.raw.get("increments", []))

    def n_nodes(self) -> int:
        return len(self.raw.get("mesh", {}).get("nodes", {}))

    def converged_u(self, increment: int = -1) -> np.ndarray:
        incs = self.increments
        if not incs:
            raise RecordError("record has no increments")
        u = incs[increment].get("u")
        if u is None:
            raise RecordError("increment %d has no converged displacement 'u'" % increment)
        return np.asarray(u, dtype=float).ravel()

    def stress_ip(self, increment: int = -1) -> Optional[Dict[int, np.ndarray]]:
        """Production per-IP stress {eid -> (8,6)} for validation, or None."""
        raw = self.increments[increment].get("stress_ip")
        if raw is None:
            return None
        return {int(eid): np.asarray(v, dtype=float) for eid, v in raw.items()}

    # -- build the neutral model ----------------------------------------- #
    def build_model(self) -> Model:
        mesh = self.raw.get("mesh", {})
        nodes = {int(k): tuple(float(c) for c in v) for k, v in mesh.get("nodes", {}).items()}
        elements: Dict[int, Element] = {}
        for k, e in mesh.get("elements", {}).items():
            eid = int(k)
            elements[eid] = Element(eid, str(e.get("type", "C3D8")),
                                    [int(n) for n in e["connectivity"]])
        model = Model(nodes=nodes, elements=elements)
        model.element_formulation = {eid: "solid_c3d8_small_strain" for eid in elements}
        model.element_material = {eid: "MAT" for eid in elements}
        model.boundaries = [self._bc(b) for b in self.raw.get("boundaries", [])]
        # named node sets, if any (lets BC targets be set names)
        model.node_sets = {str(k): [int(n) for n in v]
                           for k, v in mesh.get("node_sets", {}).items()}
        return model

    @staticmethod
    def _bc(b: Mapping[str, Any]) -> _BC:
        kind = str(b.get("kind", "value"))
        if kind.upper() in ("ENCASTRE", "PINNED", "XSYMM", "YSYMM", "ZSYMM",
                             "XASYMM", "YASYMM", "ZASYMM"):
            return _BC(target=b["target"], kind=kind.upper())
        comp = int(b.get("dof", 0))
        return _BC(target=b["target"], dof_start=comp, dof_end=comp,
                   value=float(b.get("value", 0.0)), kind="value")

    def external_load(self, dof_manager) -> np.ndarray:
        """Assemble the concentrated external force vector from record loads."""
        f = np.zeros(dof_manager.ndof, dtype=float)
        for c in self.raw.get("loads", {}).get("cload", []):
            nid, comp, val = int(c["node"]), int(c["dof"]), float(c["value"])
            f[dof_manager.node_dofs(nid)[comp - 1]] += val
        return f


# --------------------------------------------------------------------------- #
def preflight_record(record: ReplayRecord, package) -> List[str]:
    """Return a list of human-readable problems: exactly which mathematically
    required fields the record is missing for THIS material package. Empty list
    means the record is sufficient. Never raises -- it is a report.
    """
    problems: List[str] = []
    raw = record.raw

    if raw.get("schema") != SCHEMA:
        problems.append("schema must be %r, got %r" % (SCHEMA, raw.get("schema")))

    prov = record.provenance
    for k in ("model_id", "regular_hash"):
        if not prov.get(k):
            problems.append("provenance.%s is required (matched-twin check)" % k)

    if not record.kinematics:
        problems.append("kinematics is required")
    elif package is not None and record.kinematics != package.raw.get("kinematics"):
        problems.append("record kinematics %r != package kinematics %r"
                        % (record.kinematics, package.raw.get("kinematics")))

    # material.props -- must be a list of numbers (do NOT use the float() accessor,
    # which would raise on a malformed record; this is a report, not a gate)
    props = (raw.get("material", {}) or {}).get("props")
    if not isinstance(props, (list, tuple)) or not props:
        problems.append("material.props (the parameter values used in Stage 1) is required")
    elif not all(isinstance(v, (int, float)) for v in props):
        problems.append("material.props must be numbers")
    elif package is not None and len(props) != package.nprops:
        problems.append("material.props has %d values but package.nprops=%d"
                        % (len(props), package.nprops))

    mesh = raw.get("mesh", {}) or {}
    nodes = mesh.get("nodes", {}) or {}
    elements = mesh.get("elements", {}) or {}
    if not nodes:
        problems.append("mesh.nodes is required")
    if not elements:
        problems.append("mesh.elements is required")
    else:
        node_ids = set()
        for k in nodes:
            try:
                node_ids.add(int(k))
            except (TypeError, ValueError):
                problems.append("mesh.nodes has a non-integer node id %r" % (k,))
        for k, e in elements.items():
            if not isinstance(e, dict):
                problems.append("element %s is not an object" % (k,))
                continue
            if str(e.get("type", "C3D8")).upper() != "C3D8":
                problems.append("element %s: only C3D8 is supported in v1 (got %r)"
                                % (k, e.get("type")))
            conn = e.get("connectivity") or []
            try:
                missing = [n for n in conn if int(n) not in node_ids]
            except (TypeError, ValueError):
                problems.append("element %s has a non-integer connectivity entry" % (k,))
                continue
            if missing:
                problems.append("element %s references undefined nodes %r" % (k, missing[:6]))

    if not raw.get("boundaries"):
        problems.append("boundaries are required (K has no free partition without them)")

    incs = raw.get("increments") or []
    if not isinstance(incs, list) or not incs:
        problems.append("at least one converged increment is required")
        incs = []
    else:
        ndof = 3 * len(nodes)
        for i, inc in enumerate(incs):
            u = inc.get("u") if isinstance(inc, dict) else None
            if not isinstance(u, (list, tuple)):
                problems.append("increment %d: converged displacement 'u' (a list) is required" % i)
            elif ndof and len(u) != ndof:
                problems.append("increment %d: 'u' length %d != 3*nnodes=%d"
                                % (i, len(u), ndof))

    # Path-dependence: the final state alone is not enough. A single-increment
    # analysis is still valid (replay from the initial state through that one
    # increment); what is required is per-increment kinematics and the initial
    # state -- NOT two or more increments.
    nstatev = getattr(package, "nstatev", 0) if package is not None else 0
    finite = record.kinematics == "finite_strain"
    if nstatev > 0 or finite:
        for i, inc in enumerate(incs):
            if not (isinstance(inc, dict) and "kinematics_ip" in inc):
                problems.append("increment %d: per-IP 'kinematics_ip' is required for "
                                "path-dependent replay (the final state alone is "
                                "insufficient)" % i)
        if nstatev > 0 and "initial_state_ip" not in raw:
            problems.append("path-dependent replay needs 'initial_state_ip' to start from")

    return problems
