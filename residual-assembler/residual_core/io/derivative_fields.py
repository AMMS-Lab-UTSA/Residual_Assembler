"""Canonical derivative-field export: schema + loader + IP normalizer (M2).

This is the solver-neutral input the field-driven residual-method run consumes.
It is the counterpart to ``abaqus_odb_export`` (which *produces* the raw per-IP
fields from an ODB): here we *load and validate* a self-describing JSON file and
normalize its per-integration-point state into the dense ``(8, nSDV)`` arrays the
M1 engine (``core.field_sensitivity``) expects.

Canonical schema ``resasm_derivative_fields_v1``::

    {
      "metadata": {
        "schema": "resasm_derivative_fields_v1",
        "element_type": "C3D8",
        "kinematics": "small_strain",
        "voigt_order": ["11","22","33","12","13","23"],
        "integration_point_order": "abaqus_label_ascending",
        "sdv_indexing": "abaqus_1_based",
        "sdv_layout": {"ddsdde": [1,36], "parameters": {"E": [37,42], "nu": [43,48]}}
      },
      "displacements": {"<node>": [u1,u2,u3], ...},   # converged nodal U
      "statev": {"<eid>": {"<ip>": [SDV...], ...}, ...} # element -> IP -> SDV vector
    }

The ``statev`` block may also be given as a flat list per element
(``{"<eid>": [[SDV...] x 8]}``); both forms normalize to ``(8, nSDV)`` with the
integration points ordered by ascending label. Missing, duplicate, or extra
integration points are hard errors -- never silently padded.

The metadata is the SINGLE SOURCE OF TRUTH for the SDV layout: there are no
hard-coded SDV slices anywhere on this path.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

SCHEMA_ID = "resasm_derivative_fields_v1"
CANON_VOIGT_ORDER = ["11", "22", "33", "12", "13", "23"]
CANON_IP_ORDER = "abaqus_label_ascending"
CANON_SDV_INDEXING = "abaqus_1_based"
SUPPORTED_ELEMENT_TYPES = ("C3D8",)
SUPPORTED_KINEMATICS = ("small_strain",)
SUPPORTED_INTEGRATION = ("full", "selective_reduced")
N_IP = 8
_IP_LABELS = list(range(1, N_IP + 1))


class DerivativeFieldError(ValueError):
    """A user-facing problem with a derivative-field file (schema, conventions,
    or integration-point structure). Raised instead of returning a silently
    wrong / partially-padded field."""


@dataclass
class DerivativeFields:
    """A validated, normalized derivative-field export.

    Attributes
    ----------
    metadata : the file's ``metadata`` block (schema, conventions, sdv_layout).
    statev : ``{eid(int) -> (8, nSDV) ndarray}`` per-IP state, IPs 1..8 ascending.
    displacements : ``{node(int) -> (3,) ndarray}`` converged nodal displacement.
    """
    metadata: Dict[str, Any]
    statev: Dict[int, np.ndarray]
    displacements: Dict[int, np.ndarray] = field(default_factory=dict)

    @property
    def sdv_layout(self) -> Dict[str, Any]:
        return self.metadata["sdv_layout"]

    @property
    def parameters(self) -> List[str]:
        return list(self.metadata.get("sdv_layout", {}).get("parameters", {}).keys())

    @property
    def element_type(self) -> str:
        return self.metadata.get("element_type", "")

    @property
    def kinematics(self) -> str:
        return self.metadata.get("kinematics", "")

    @property
    def integration(self) -> str:
        """Volumetric integration scheme: 'full' (default) or 'selective_reduced'
        (Abaqus C3D8 mean-dilatation B-bar)."""
        return str(self.metadata.get("integration", "full"))


# --------------------------------------------------------------------------- #
# integration-point normalization
# --------------------------------------------------------------------------- #
def _normalize_element_statev(eid, per_ip) -> np.ndarray:
    """One element's per-IP SDV vectors -> (8, nSDV), IPs ordered 1..8 ascending.

    Accepts the canonical nested mapping ``{ip_label -> [SDV...]}`` or a flat
    list ``[[SDV...] x 8]`` (already in ascending IP order). Rejects missing,
    duplicate, or extra integration points and ragged SDV rows.
    """
    if isinstance(per_ip, dict):
        labels: List[int] = []
        for ip_key in per_ip:
            try:
                labels.append(int(ip_key))
            except (TypeError, ValueError):
                raise DerivativeFieldError(
                    "element %r has a non-integer integration-point label %r"
                    % (eid, ip_key))
        # a JSON object cannot carry duplicate keys, but a normalized-from-ints
        # form could; guard anyway.
        if len(set(labels)) != len(labels):
            raise DerivativeFieldError(
                "element %r has duplicate integration-point labels %r" % (eid, labels))
        ordered = sorted(labels)
        if ordered != _IP_LABELS:
            raise DerivativeFieldError(
                "element %r integration-point labels must be exactly 1..8 ascending, "
                "got %r (missing %r, extra %r)"
                % (eid, ordered, sorted(set(_IP_LABELS) - set(ordered)),
                   sorted(set(ordered) - set(_IP_LABELS))))
        rows = [np.asarray(per_ip[k], dtype=float)
                for k in _keys_sorted_by_int(per_ip)]
    elif isinstance(per_ip, (list, tuple)):
        if len(per_ip) != N_IP:
            raise DerivativeFieldError(
                "element %r must have exactly %d integration points, got %d"
                % (eid, N_IP, len(per_ip)))
        rows = [np.asarray(v, dtype=float) for v in per_ip]
    else:
        raise DerivativeFieldError(
            "element %r statev must be a mapping {ip: [SDV...]} or a list of 8 "
            "SDV vectors, got %s" % (eid, type(per_ip).__name__))

    nsdv = rows[0].shape[0] if rows[0].ndim == 1 else -1
    for ip, r in enumerate(rows, start=1):
        if r.ndim != 1:
            raise DerivativeFieldError(
                "element %r IP %d SDV entry must be a flat vector, got shape %r"
                % (eid, ip, tuple(r.shape)))
        if r.shape[0] != nsdv:
            raise DerivativeFieldError(
                "element %r has ragged SDV vectors: IP1 has %d values but IP %d "
                "has %d" % (eid, nsdv, ip, r.shape[0]))
    return np.array(rows, dtype=float)


def _keys_sorted_by_int(d):
    return sorted(d.keys(), key=lambda k: int(k))


def normalize_statev(raw: Dict[Any, Any]) -> Dict[int, np.ndarray]:
    """Normalize a raw ``statev`` block to ``{eid(int) -> (8, nSDV) ndarray}``."""
    if not isinstance(raw, dict) or not raw:
        raise DerivativeFieldError(
            "'statev' must be a non-empty mapping {element_id -> per-IP SDV data}")
    out: Dict[int, np.ndarray] = {}
    for eid_key, per_ip in raw.items():
        try:
            eid = int(eid_key)
        except (TypeError, ValueError):
            raise DerivativeFieldError(
                "element id %r is not an integer" % (eid_key,))
        out[eid] = _normalize_element_statev(eid, per_ip)
    return out


# --------------------------------------------------------------------------- #
# metadata validation
# --------------------------------------------------------------------------- #
def _validate_metadata(meta: Dict[str, Any]) -> None:
    if not isinstance(meta, dict):
        raise DerivativeFieldError("'metadata' block is missing or not a mapping")
    schema = meta.get("schema")
    if schema != SCHEMA_ID:
        raise DerivativeFieldError(
            "unsupported derivative-field schema %r; expected %r" % (schema, SCHEMA_ID))

    et = str(meta.get("element_type", "")).upper()
    if et not in SUPPORTED_ELEMENT_TYPES:
        raise DerivativeFieldError(
            "unsupported element_type %r; the field-driven path supports %s only"
            % (meta.get("element_type"), "/".join(SUPPORTED_ELEMENT_TYPES)))

    kin = str(meta.get("kinematics", "")).lower()
    if kin not in SUPPORTED_KINEMATICS:
        raise DerivativeFieldError(
            "unsupported kinematics %r; only %s is supported (finite-strain is a "
            "later milestone)" % (meta.get("kinematics"), "/".join(SUPPORTED_KINEMATICS)))

    integ = str(meta.get("integration", "full")).lower()
    if integ not in SUPPORTED_INTEGRATION:
        raise DerivativeFieldError(
            "unsupported integration %r; use one of %s (Abaqus C3D8 -> "
            "'selective_reduced')" % (meta.get("integration"),
                                      "/".join(SUPPORTED_INTEGRATION)))

    voigt = meta.get("voigt_order")
    if list(voigt or []) != CANON_VOIGT_ORDER:
        raise DerivativeFieldError(
            "voigt_order must be %r (Abaqus order), got %r" % (CANON_VOIGT_ORDER, voigt))

    ip_order = str(meta.get("integration_point_order", ""))
    if ip_order != CANON_IP_ORDER:
        raise DerivativeFieldError(
            "integration_point_order must be %r, got %r" % (CANON_IP_ORDER, ip_order))

    sdv_indexing = str(meta.get("sdv_indexing", ""))
    if sdv_indexing != CANON_SDV_INDEXING:
        raise DerivativeFieldError(
            "sdv_indexing must be %r (1-based inclusive Abaqus SDV numbers), got %r"
            % (CANON_SDV_INDEXING, sdv_indexing))

    if not isinstance(meta.get("sdv_layout"), dict):
        raise DerivativeFieldError("metadata is missing the 'sdv_layout' block")
    if "ddsdde" not in meta["sdv_layout"]:
        raise DerivativeFieldError("sdv_layout is missing the 'ddsdde' range")


# --------------------------------------------------------------------------- #
# loading
# --------------------------------------------------------------------------- #
def load_derivative_fields(source: Any) -> DerivativeFields:
    """Load and validate a derivative-field file (path) or an already-parsed dict.

    Raises :class:`DerivativeFieldError` on any schema, convention, or
    integration-point problem.
    """
    if isinstance(source, str):
        try:
            with open(source, "r", encoding="utf-8") as fh:
                data = json.load(fh)
        except (OSError, json.JSONDecodeError) as exc:
            raise DerivativeFieldError("cannot read derivative-field file %r: %s"
                                       % (source, exc))
    elif isinstance(source, dict):
        data = source
    else:
        raise DerivativeFieldError(
            "load_derivative_fields expects a path or a dict, got %s"
            % type(source).__name__)

    if not isinstance(data, dict):
        raise DerivativeFieldError("derivative-field file must be a JSON object")

    _validate_metadata(data.get("metadata"))
    if "statev" not in data:
        raise DerivativeFieldError("derivative-field file is missing the 'statev' block")

    statev = normalize_statev(data["statev"])

    disp_raw = data.get("displacements") or {}
    displacements: Dict[int, np.ndarray] = {}
    for node_key, vec in disp_raw.items():
        try:
            nid = int(node_key)
        except (TypeError, ValueError):
            raise DerivativeFieldError("displacement node id %r is not an integer"
                                       % (node_key,))
        displacements[nid] = np.asarray(vec, dtype=float).ravel()

    return DerivativeFields(metadata=data["metadata"], statev=statev,
                            displacements=displacements)


def is_derivative_field_dict(data: Any) -> bool:
    """True if ``data`` looks like a resasm_derivative_fields_v1 export."""
    return (isinstance(data, dict) and isinstance(data.get("metadata"), dict)
            and data["metadata"].get("schema") == SCHEMA_ID)


def displacements_to_vector(displacements: Dict[int, np.ndarray],
                            dof_manager) -> np.ndarray:
    """Scatter ``{node -> (3,)}`` onto a global (ndof,) vector in the DofManager's
    node-major (UX,UY,UZ) order. Nodes without a supplied displacement stay zero
    (the small-strain field-driven K and R_,a do not depend on U, so this only
    affects provenance / future output chain-rules)."""
    u = np.zeros(dof_manager.ndof, dtype=float)
    for nid, vec in displacements.items():
        if nid not in dof_manager.node_index:
            continue
        gdofs = dof_manager.node_dofs(nid)
        for k, val in enumerate(vec):
            if k < len(gdofs):
                u[gdofs[k]] = float(val)
    return u
