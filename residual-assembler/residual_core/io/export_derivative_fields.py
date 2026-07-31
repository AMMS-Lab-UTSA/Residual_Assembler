#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""export_derivative_fields.py -- ODB -> resasm_derivative_fields_v1 (M3).

Turns a real Abaqus analysis into the canonical derivative-field JSON that the
M2 recipe (`resasm run sensitivity.yaml`) already consumes. It reuses the good
per-integration-point extraction in ``residual_core/io/abaqus_odb_export.py``
(NOT the flat ``scripts/extract_odb_fields.py``).

Run under Abaqus' Python (Python 2.7 + odbAccess)::

    abaqus python export_derivative_fields.py -- \
        --odb nonuniform_c3d8.odb --step Step-1 --frame -1 \
        --layout derivative_layout.json --output derivative_fields.json

The SDV semantics (which SDV slots are DDSDDE, which are d sigma/d a_i) are NOT
inferable from an ODB. They come from an external ``resasm_sdv_layout_v1`` sidecar
(``--layout``); that layout is the single source of truth and is copied verbatim
into the output's metadata.

This module is import-safe under normal CPython 3 too: the odbAccess-dependent
code lives in ``main``/``harvest_*``; the pure builder/validator
(``build_derivative_fields``, ``load_sdv_layout``) has no Abaqus dependency and is
unit-tested offline.

Precision note: Abaqus stores field output (including SDV and U) in the ODB as
SINGLE precision. Values round-tripped through the ODB are therefore accurate to
~1e-7 relative, not to double precision -- this bounds every downstream check.
"""
from __future__ import print_function

import json
import os
import sys

SCHEMA_DERIVATIVE_FIELDS = "resasm_derivative_fields_v1"
SCHEMA_SDV_LAYOUT = "resasm_sdv_layout_v1"
N_IP = 8
_IP_LABELS = list(range(1, N_IP + 1))
_LAYOUT_META_KEYS = ("element_type", "kinematics", "voigt_order",
                     "integration_point_order", "sdv_indexing", "sdv_layout")


class ExportError(Exception):
    """A user-facing problem exporting a derivative-field file."""


# --------------------------------------------------------------------------- #
# sidecar SDV layout (resasm_sdv_layout_v1)
# --------------------------------------------------------------------------- #
def load_sdv_layout(source):
    """Load + validate a ``resasm_sdv_layout_v1`` sidecar (path or dict).

    Returns a metadata dict carrying element_type, kinematics, voigt_order,
    integration_point_order, sdv_indexing and sdv_layout -- ready to embed in the
    derivative-fields metadata. Raises :class:`ExportError` on a bad schema.
    """
    if isinstance(source, str):
        with open(source, "r") as fh:
            data = json.load(fh)
    elif isinstance(source, dict):
        data = source
    else:
        raise ExportError("layout must be a path or dict, got %s"
                          % type(source).__name__)
    if data.get("schema") != SCHEMA_SDV_LAYOUT:
        raise ExportError("layout schema must be %r, got %r"
                          % (SCHEMA_SDV_LAYOUT, data.get("schema")))
    for k in _LAYOUT_META_KEYS:
        if k not in data:
            raise ExportError("layout is missing required key %r" % k)
    layout = data["sdv_layout"]
    if not isinstance(layout, dict) or "ddsdde" not in layout:
        raise ExportError("layout.sdv_layout must be a mapping with a 'ddsdde' range")
    meta = dict((k, data[k]) for k in _LAYOUT_META_KEYS)
    # optional: volumetric integration scheme (Abaqus C3D8 -> selective_reduced)
    if "integration" in data:
        meta["integration"] = data["integration"]
    return meta


def _max_sdv_needed(sdv_layout):
    """Largest 1-based SDV number the layout references."""
    hi = int(sdv_layout["ddsdde"][1])
    for rng in dict(sdv_layout.get("parameters", {})).values():
        hi = max(hi, int(rng[1]))
    return hi


# --------------------------------------------------------------------------- #
# pure builder (no Abaqus) -- unit-tested offline
# --------------------------------------------------------------------------- #
def build_derivative_fields(displacements, statev_nested, layout_meta,
                            provenance=None, mesh_element_ids=None):
    """Assemble the ``resasm_derivative_fields_v1`` document (a plain dict).

    Parameters
    ----------
    displacements : ``{node_label(str/int) -> [u1,u2,u3]}``.
    statev_nested : ``{eid(str/int) -> {ip_label(int) -> [SDV...]}}`` dense per-IP
        SDV vectors (index i holds SDV(i+1)).
    layout_meta : the dict from :func:`load_sdv_layout`.
    provenance : optional dict merged into metadata (odb, step, frame, times).
    mesh_element_ids : optional iterable of every element id that must be present
        (so a missing element is a hard error, not a silent gap).

    Raises :class:`ExportError` on any integration-point / coverage / SDV-length
    problem.
    """
    max_sdv = _max_sdv_needed(layout_meta["sdv_layout"])

    statev_out = {}
    for eid, per_ip in statev_nested.items():
        if not isinstance(per_ip, dict):
            raise ExportError("element %r statev must be an {ip: vector} mapping" % (eid,))
        labels = sorted(int(k) for k in per_ip)
        if labels != _IP_LABELS:
            raise ExportError(
                "element %r integration-point labels must be exactly 1..8, got %r"
                % (eid, labels))
        rows = {}
        length = None
        for ip in _IP_LABELS:
            vec = list(per_ip[ip]) if ip in per_ip else list(per_ip[str(ip)])
            if length is None:
                length = len(vec)
            elif len(vec) != length:
                raise ExportError(
                    "element %r has ragged SDV vectors (IP1 len %d, IP %d len %d)"
                    % (eid, length, ip, len(vec)))
            rows[str(ip)] = [float(x) for x in vec]
        if length < max_sdv:
            raise ExportError(
                "element %r exports %d SDVs but the layout needs at least %d "
                "(is *Depvar large enough and are all SDVk requested in "
                "*Element Output?)" % (eid, length, max_sdv))
        statev_out[str(int(eid))] = rows

    if mesh_element_ids is not None:
        want = set(str(int(e)) for e in mesh_element_ids)
        got = set(statev_out.keys())
        missing = sorted(want - got)
        if missing:
            raise ExportError("no exported SDV field for element(s) %r present in "
                              "the mesh" % missing[:8])

    disp_out = {}
    for node, vec in displacements.items():
        disp_out[str(int(node))] = [float(x) for x in vec]

    metadata = {"schema": SCHEMA_DERIVATIVE_FIELDS}
    for k in _LAYOUT_META_KEYS:
        metadata[k] = layout_meta[k]
    if "integration" in layout_meta:
        metadata["integration"] = layout_meta["integration"]
    if provenance:
        metadata["provenance"] = dict(provenance)

    return {"metadata": metadata, "displacements": disp_out, "statev": statev_out}


# --------------------------------------------------------------------------- #
# Abaqus harvest (odbAccess) -- only runs under Abaqus python
# --------------------------------------------------------------------------- #
def _import_odb_helpers():
    """Import the shared ODB extraction helpers as a standalone module (they are
    py2.7-safe and dependency-free, so this works under Abaqus python without
    importing the whole residual_core package)."""
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    import abaqus_odb_export as odbx           # noqa: E402
    return odbx


def harvest_statev_nested(frame, region, instance):
    """``{eid(str) -> {ip(int) -> [SDV1..SDVmax]}}`` with IP labels preserved.

    Dense per-IP vector: index i holds SDV(i+1); gaps (absent SDVk) are 0.0.
    """
    from abaqusConstants import INTEGRATION_POINT      # Abaqus-only
    keys = [k for k in frame.fieldOutputs.keys()
            if k.upper().startswith("SDV") and k[3:].isdigit()]
    if not keys:
        raise ExportError("the ODB has no SDV field outputs; request SDV in "
                          "*Element Output and ensure *Depvar is set")
    nums = sorted(int(k[3:]) for k in keys)
    maxn = nums[-1]
    reg = region if region is not None else instance

    acc = {}
    for n in nums:
        fld = frame.fieldOutputs["SDV%d" % n].getSubset(
            region=reg, position=INTEGRATION_POINT)
        for v in fld.values:
            el = str(int(v.elementLabel))
            ip = int(v.integrationPoint)
            row = acc.setdefault(el, {}).setdefault(ip, [0.0] * maxn)
            row[n - 1] = float(v.data)
    return acc


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        from odbAccess import openOdb
    except Exception as exc:                            # not under Abaqus python
        sys.stderr.write("ERROR: run with Abaqus python (odbAccess): %s\n" % exc)
        return 2

    opts = _parse_args(argv)
    odbx = _import_odb_helpers()

    layout_meta = load_sdv_layout(opts["layout"])

    odb = openOdb(opts["odb"], readOnly=True)
    try:
        step_name, step = odbx.pick_step(odb, opts["step"])
        frame_idx = _resolve_frame(step, opts["frame"])
        frame = step.frames[frame_idx]
        inst_name, instance = odbx.pick_instance(odb.rootAssembly, opts["instance"])
        region = odbx.elset_region(instance, odb, opts["elset"])

        displacements = odbx.nodal_field(frame, "U", instance)
        if displacements is None:
            raise ExportError("no nodal displacement field 'U' in the selected frame")
        statev_nested = harvest_statev_nested(frame, region, instance)
        _nodes, mesh_elems = odbx.collect_mesh(instance, region)

        provenance = {
            "odb": os.path.basename(opts["odb"]),
            "step": step_name,
            "frame": frame_idx,
            "step_time": float(frame.frameValue),
            "total_time": float(getattr(step, "totalTime", 0.0)) + float(frame.frameValue),
        }
        doc = build_derivative_fields(displacements, statev_nested, layout_meta,
                                      provenance=provenance,
                                      mesh_element_ids=mesh_elems.keys())
    finally:
        odb.close()

    with open(opts["output"], "w") as fh:
        json.dump(doc, fh, indent=2)
    print("wrote %s (%d elements, %d nodes, step=%s frame=%d)"
          % (opts["output"], len(doc["statev"]), len(doc["displacements"]),
             provenance["step"], provenance["frame"]))
    return 0


def _resolve_frame(step, spec):
    nframes = len(step.frames)
    if nframes == 0:
        raise ExportError("selected step has no frames")
    j = int(spec)
    if j < 0:
        j += nframes
    if j < 0 or j >= nframes:
        raise ExportError("frame index %s out of range 0..%d" % (spec, nframes - 1))
    return j


def _parse_args(argv):
    if "--" in argv:
        argv = argv[argv.index("--") + 1:]
    opts = {"odb": None, "step": None, "frame": "-1", "instance": None,
            "elset": None, "layout": None, "output": "derivative_fields.json"}
    i = 0
    while i < len(argv):
        a = argv[i]
        key = a[2:] if a.startswith("--") else None
        if key is None or key not in opts:
            sys.stderr.write("ERROR: unknown/misplaced argument: %s\n" % a)
            sys.exit(2)
        if i + 1 >= len(argv):
            sys.stderr.write("ERROR: missing value for --%s\n" % key)
            sys.exit(2)
        opts[key] = argv[i + 1]
        i += 2
    for req in ("odb", "layout"):
        if not opts[req]:
            sys.stderr.write("ERROR: --%s is required\n" % req)
            sys.exit(2)
    return opts


if __name__ == "__main__":
    sys.exit(main())
