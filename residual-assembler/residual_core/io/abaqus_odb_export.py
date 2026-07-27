#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extract_abaqus_fields.py -- Abaqus/Python ODB field extractor.

Writes the `fields.json` exchange format defined in residual_core/CONTRACT.md
section 5 (undeformed nodal coords, C3D8 connectivity, and per-frame nodal U,
nodal RF, integration-point S, and optional SDV) so the external residual
assembler and the UMAT replay adapter can consume Abaqus results.

RUN IT WITH ABAQUS' PYTHON (this is Python 2.7 + odbAccess; it will NOT run
under a normal python3).  Examples:

    abaqus python extract_abaqus_fields.py -- \
        --odb Compression111.odb --out fields.json

    abaqus python extract_abaqus_fields.py -- \
        --odb Job-1.odb --instance PART-1-1 --step Step-1 \
        --frames last --elset Set-1 --out fields.json

    # all frames, comma list, or explicit indices:
    abaqus python extract_abaqus_fields.py -- --odb job.odb --frames all
    abaqus python extract_abaqus_fields.py -- --odb job.odb --frames 0,5,10

    # via CAE (no GUI):
    abaqus cae noGUI=extract_abaqus_fields.py -- --odb job.odb --out fields.json

CLI:
    --odb       (required) path to the .odb
    --instance  instance name (default: first instance in the assembly)
    --step      step name    (default: last step)
    --frames    'last' (default) | 'all' | comma list of frame indices
    --elset     element-set name to restrict S/SDV/elements (optional)
    --out       output json path (default: fields.json)

Notes
-----
* S is read at position=INTEGRATION_POINT and grouped by element label, then
  sorted by integrationPoint so index 0..7 == Abaqus IP 1..8 (CONTRACT sec 2).
* S component order is Abaqus order (S11,S22,S33,S12,S13,S23) -- that is the
  native ODB order for C3D8.
* SDVs appear in the ODB as separate scalar fields 'SDV1','SDV2',...  This
  script exports whichever SDVk exist and records their numbers in the extra
  top-level key 'sdv_labels' (a CONTRACT-compatible additive extension used by
  umat_replay.py to line SDVs up with STATEV indices).
"""
from __future__ import print_function

import sys
import os
import json

try:
    from odbAccess import openOdb
    from abaqusConstants import INTEGRATION_POINT
except Exception as exc:  # not running under Abaqus python
    sys.stderr.write(
        "ERROR: this script must be run with Abaqus' Python (odbAccess).\n"
        "       e.g.  abaqus python extract_abaqus_fields.py -- --odb job.odb\n"
        "       import error: %s\n" % exc)
    sys.exit(2)


def _f(x):
    """Coerce an ODB scalar (possibly numpy) to a plain python float."""
    return float(x)


def _vec(data):
    """Coerce an ODB .data tuple/array to a list of plain floats."""
    return [float(v) for v in data]


# --------------------------------------------------------------------------
# tiny argparse-free CLI (argparse exists in 2.7 but keep this dependency-free
# and robust to the Abaqus '--' separator)
# --------------------------------------------------------------------------
def parse_args(argv):
    # drop everything up to and including a lone '--' if Abaqus inserted one
    if "--" in argv:
        argv = argv[argv.index("--") + 1:]
    opts = {
        "odb": None, "instance": None, "step": None,
        "frames": "last", "elset": None, "out": "fields.json",
    }
    i = 0
    while i < len(argv):
        a = argv[i]
        key = a[2:] if a.startswith("--") else None
        if key is None or key not in opts:
            sys.stderr.write("ERROR: unknown/misplaced argument: %s\n" % a)
            sys.stderr.write(__doc__)
            sys.exit(2)
        if i + 1 >= len(argv):
            sys.stderr.write("ERROR: missing value for --%s\n" % key)
            sys.exit(2)
        opts[key] = argv[i + 1]
        i += 2
    if not opts["odb"]:
        sys.stderr.write("ERROR: --odb is required.\n")
        sys.stderr.write(__doc__)
        sys.exit(2)
    return opts


def pick_instance(assembly, name):
    inst_names = list(assembly.instances.keys())
    if not inst_names:
        raise RuntimeError("no instances found in the assembly.")
    if name is None:
        chosen = inst_names[0]
        if len(inst_names) > 1:
            sys.stderr.write("WARNING: multiple instances %s; using '%s'. "
                             "Use --instance to choose.\n" % (inst_names, chosen))
        return chosen, assembly.instances[chosen]
    if name not in assembly.instances:
        raise RuntimeError("instance '%s' not found. Available: %s"
                           % (name, inst_names))
    return name, assembly.instances[name]


def pick_step(odb, name):
    step_names = list(odb.steps.keys())
    if not step_names:
        raise RuntimeError("no steps found in the odb.")
    if name is None:
        chosen = step_names[-1]
        return chosen, odb.steps[chosen]
    if name not in odb.steps:
        raise RuntimeError("step '%s' not found. Available: %s"
                           % (name, step_names))
    return name, odb.steps[name]


def pick_frames(step, spec):
    nframes = len(step.frames)
    if nframes == 0:
        raise RuntimeError("step has no frames.")
    spec = (spec or "last").strip().lower()
    if spec == "last":
        idxs = [nframes - 1]
    elif spec == "all":
        idxs = list(range(nframes))
    else:
        idxs = []
        for tok in spec.split(","):
            tok = tok.strip()
            if tok == "":
                continue
            j = int(tok)
            if j < 0:
                j += nframes
            if j < 0 or j >= nframes:
                raise RuntimeError("frame index %s out of range 0..%d"
                                   % (tok, nframes - 1))
            idxs.append(j)
    return idxs


def elset_region(instance, odb, elset):
    """Return an OdbSet region for `elset`, or None (whole instance)."""
    if not elset:
        return None
    if elset in instance.elementSets:
        return instance.elementSets[elset]
    if elset in odb.rootAssembly.elementSets:
        return odb.rootAssembly.elementSets[elset]
    raise RuntimeError("elset '%s' not found in instance or assembly." % elset)


def collect_mesh(instance, region):
    """Undeformed coords + C3D8 connectivity (optionally restricted to region)."""
    if region is not None and hasattr(region, "elements") and region.elements:
        elems = region.elements
        # region.elements may be nested per-instance
        if elems and not hasattr(elems[0], "label"):
            elems = elems[0]
    else:
        elems = instance.elements

    elements = {}
    used_nodes = set()
    skipped = 0
    for e in elems:
        etype = str(e.type)
        if not etype.startswith("C3D8"):
            skipped += 1
            continue
        conn = [int(n) for n in e.connectivity]
        elements[str(int(e.label))] = conn
        for n in conn:
            used_nodes.add(int(n))
    if skipped:
        sys.stderr.write("WARNING: skipped %d non-C3D8 elements.\n" % skipped)

    nodes = {}
    restrict = len(used_nodes) > 0 and region is not None
    for nd in instance.nodes:
        lab = int(nd.label)
        if restrict and lab not in used_nodes:
            continue
        nodes[str(lab)] = [_f(c) for c in nd.coordinates]
    return nodes, elements


def nodal_field(frame, key, instance):
    """Return {nodeLabel: [x,y,z]} for a nodal field, or None if absent."""
    if key not in frame.fieldOutputs:
        return None
    fld = frame.fieldOutputs[key].getSubset(region=instance)
    out = {}
    for v in fld.values:
        out[str(int(v.nodeLabel))] = _vec(v.data)
    return out


def ip_stress(frame, region, instance):
    """Return {elementLabel: [[6 comps] x 8 IPs]} for S, IPs sorted ascending."""
    if "S" not in frame.fieldOutputs:
        return None
    fld = frame.fieldOutputs["S"]
    # restrict to the elset region if given, else to the picked instance so a
    # multi-instance ODB does not merge different instances' element label "1"
    # (matches the nodal U/RF restriction in nodal_field).
    reg = region if region is not None else instance
    fld = fld.getSubset(region=reg, position=INTEGRATION_POINT)
    tmp = {}
    for v in fld.values:
        el = str(int(v.elementLabel))
        tmp.setdefault(el, []).append((int(v.integrationPoint), _vec(v.data)))
    out = {}
    for el, pairs in tmp.items():
        pairs.sort(key=lambda p: p[0])          # IP 1..8 ascending
        out[el] = [d for (_ip, d) in pairs]
    return out


def sdv_fields(frame, region, instance):
    """
    Return ({elementLabel: [[k SDVs] x 8 IPs]}, sdv_labels) or (None, []).
    SDVs are separate scalar fields 'SDV1'.. in the ODB.
    """
    keys = [k for k in frame.fieldOutputs.keys()
            if k.upper().startswith("SDV") and k[3:].isdigit()]
    if not keys:
        return None, []
    keys.sort(key=lambda k: int(k[3:]))
    labels = [int(k[3:]) for k in keys]

    # restrict to the elset region if given, else to the picked instance so a
    # multi-instance ODB does not merge different instances' element label "1".
    reg = region if region is not None else instance

    # per element -> per IP index -> list of SDV values in `keys` order
    acc = {}
    for col, key in enumerate(keys):
        fld = frame.fieldOutputs[key].getSubset(
            region=reg, position=INTEGRATION_POINT)
        for v in fld.values:
            el = str(int(v.elementLabel))
            ip = int(v.integrationPoint)
            acc.setdefault(el, {}).setdefault(ip, {})[col] = _f(v.data)

    out = {}
    ncol = len(keys)
    for el, byip in acc.items():
        rows = []
        for ip in sorted(byip.keys()):
            row = byip[ip]
            rows.append([row.get(c, 0.0) for c in range(ncol)])
        out[el] = rows
    return out, labels


def main(argv):
    opts = parse_args(argv)
    if not os.path.exists(opts["odb"]):
        sys.stderr.write("ERROR: odb not found: %s\n" % opts["odb"])
        sys.exit(3)

    print("opening odb: %s" % opts["odb"])
    odb = openOdb(opts["odb"], readOnly=True)
    try:
        inst_name, instance = pick_instance(odb.rootAssembly, opts["instance"])
        step_name, step = pick_step(odb, opts["step"])
        region = elset_region(instance, odb, opts["elset"])
        frame_idxs = pick_frames(step, opts["frames"])
        print("instance=%s  step=%s  frames=%s  elset=%s"
              % (inst_name, step_name, frame_idxs, opts["elset"]))

        nodes, elements = collect_mesh(instance, region)
        print("mesh: %d nodes, %d C3D8 elements" % (len(nodes), len(elements)))

        frames_out = []
        sdv_labels_seen = []
        for j in frame_idxs:
            frame = step.frames[j]
            rec = {
                "frame": int(frame.frameId) if hasattr(frame, "frameId") else j,
                "step": step_name,
                "time": _f(frame.frameValue),
            }
            U = nodal_field(frame, "U", instance)
            if U is not None:
                rec["U"] = U
            RF = nodal_field(frame, "RF", instance)
            if RF is not None:
                rec["RF"] = RF
            S = ip_stress(frame, region, instance)
            if S is not None:
                rec["S"] = S
            SDV, labels = sdv_fields(frame, region, instance)
            if SDV is not None:
                rec["SDV"] = SDV
                if labels and not sdv_labels_seen:
                    sdv_labels_seen = labels
            frames_out.append(rec)
            print("  frame %d (t=%.6g): U=%s RF=%s S=%s SDV=%s"
                  % (j, rec["time"], "U" in rec, "RF" in rec, "S" in rec, "SDV" in rec))

        out = {
            "odb": os.path.basename(opts["odb"]),
            "instance": inst_name,
            "element_type": "C3D8",
            "nodes": nodes,
            "elements": elements,
            "frames": frames_out,
        }
        if sdv_labels_seen:
            out["sdv_labels"] = sdv_labels_seen

        with open(opts["out"], "w") as fh:
            json.dump(out, fh, separators=(",", ":"), sort_keys=True)
        print("wrote %s  (%d frame record(s)%s)"
              % (opts["out"], len(frames_out),
                 ", sdv_labels=%s" % sdv_labels_seen if sdv_labels_seen else ""))
    finally:
        odb.close()


if __name__ == "__main__":
    main(sys.argv[1:])
