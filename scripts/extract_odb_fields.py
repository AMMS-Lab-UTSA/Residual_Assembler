"""extract_odb_fields.py — export ODB fields to a solver-neutral JSON.

Step 2 of the roadmap. Reads integration-point stress (S), reaction forces (RF),
displacements (U) and state variables (SDV) from an Abaqus ODB and writes the
``fields.json`` consumed by ``resasm assemble --mode stress-driven`` and by
``compare_residuals.py``.

Must be run with the **Abaqus Python** interpreter (it imports ``odbAccess``):
    abaqus python scripts/extract_odb_fields.py --odb my_job.odb --out fields.json

When ``odbAccess`` is unavailable (i.e. not inside Abaqus), it prints
"Abaqus not available: validation pending" and exits 0.

Output schema (matches attach_results / neutral_model_io field_refs):
    {
      "stress_ip": {"<eid>": [[s11,s22,s33,s12,s13,s23], ... per IP]},
      "reactions": {"<node>": [rf1, rf2, rf3]},
      "displacements": {"<node>": [u1, u2, u3]},
      "statev": {"<eid>": [[sdv...] per IP]},
      "meta": {"odb": "...", "step": "...", "frame": <int>}
    }
"""

from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from _abaqus_env import skip_if_no_abaqus, SKIP_MESSAGE  # noqa: E402


def extract(odb_path, step_name=None, frame_index=-1):
    import odbAccess  # Abaqus-only

    odb = odbAccess.openOdb(path=odb_path, readOnly=True)
    step = (odb.steps[step_name] if step_name else
            odb.steps[list(odb.steps.keys())[-1]])
    frame = step.frames[frame_index]

    out = {"stress_ip": {}, "reactions": {}, "displacements": {},
           "statev": {}, "meta": {"odb": odb_path, "step": step.name,
                                  "frame": frame.frameId}}

    fo = frame.fieldOutputs
    if "S" in fo:
        for v in fo["S"].values:
            eid = str(v.elementLabel)
            out["stress_ip"].setdefault(eid, []).append(
                [float(c) for c in v.data])   # (s11,s22,s33,s12,s13,s23)
    if "RF" in fo:
        for v in fo["RF"].values:
            out["reactions"][str(v.nodeLabel)] = [float(c) for c in v.data]
    if "U" in fo:
        for v in fo["U"].values:
            out["displacements"][str(v.nodeLabel)] = [float(c) for c in v.data]
    # any SDVn present -> collect per element/IP
    sdv_keys = sorted(k for k in fo.keys() if k.startswith("SDV"))
    if sdv_keys:
        per = {}
        for k in sdv_keys:
            for v in fo[k].values:
                per.setdefault(str(v.elementLabel), []).append(float(v.data))
        out["statev"] = per

    odb.close()
    return out


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Export ODB fields to JSON (skips w/o Abaqus)")
    p.add_argument("--odb", required=True)
    p.add_argument("--out", default="fields.json")
    p.add_argument("--step", default=None)
    p.add_argument("--frame", type=int, default=-1)
    args = p.parse_args(argv)

    # graceful skip when not inside Abaqus python
    skip_if_no_abaqus(need_odb=True)

    if not os.path.exists(args.odb):
        print("%s (odb not found: %s)" % (SKIP_MESSAGE, args.odb))
        return 0

    data = extract(args.odb, args.step, args.frame)
    with open(args.out, "w") as fh:
        json.dump(data, fh, indent=2)
    print("wrote %s (%d elements with stress, %d reactions)"
          % (args.out, len(data["stress_ip"]), len(data["reactions"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
