#!/usr/bin/env python3
"""PROGRAM 2 (collaborator) demo: sensitivities from a JHU package + a saved run.

This is the collaborator workflow. The tool consumes ONLY the packaged binary +
manifest and a saved replay record; it never sees JHU source.

For a self-contained demo we first SIMULATE what a collaborator would already
have -- JHU's shipped package (opaque .so + manifest) and a saved production
analysis (base + p+/-dp records) -- using the isolated test fixture. Then the
real work is the single ``run_request`` call and the finite-difference check
against the perturbed records.

    python examples/replay_elastic_c3d8/run_replay.py            # displacement-driven
    python examples/replay_elastic_c3d8/run_replay.py --drive force

Needs gfortran once (to compile the opaque demo binary). Writes the four contract
documents next to this file.
"""

import argparse
import json
import os
import shutil
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))
for p in (_ROOT, os.path.join(_ROOT, "tests", "fixtures", "replay_elastic")):
    if p not in sys.path:
        sys.path.insert(0, p)

from residual_core.replay import run_request          # the Program-2 application
import make_fixtures                                   # demo scaffolding only


def _vm(sig):
    s11, s22, s33, s12, s13, s23 = sig
    return float(np.sqrt(0.5 * ((s11 - s22) ** 2 + (s22 - s33) ** 2 + (s33 - s11) ** 2)
                         + 3 * (s12 * s12 + s13 * s13 + s23 * s23)))


def _q(rec, req):
    inc = rec["increments"][0]
    if req["type"] == "displacement":
        return inc["u"][3 * (req["node"] - 1) + req["dof"] - 1]
    if req["type"] == "reaction":
        return inc["reactions"][3 * (req["node"] - 1) + req["dof"] - 1]
    sig = np.array(inc["stress_ip"][str(req["element"])][req["ip"]])
    return _vm(sig) if req.get("von_mises") else sig[req["component"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drive", choices=["force", "disp"], default="disp")
    args = ap.parse_args()
    if shutil.which("gfortran") is None:
        print("  gfortran needed to build the opaque demo binary."); return 1

    # ---- simulated: what JHU shipped + the collaborator's saved analyses ----
    fx = make_fixtures.build_fixtures(_HERE, drive=args.drive)
    request = {
        "schema": "resasm_sensitivity_request_v1",
        "parameters": ["E", "nu"],
        "outputs": [
            {"type": "displacement", "node": 11, "dof": 1},
            {"type": "reaction", "node": 11, "dof": 3},
            {"type": "stress", "element": 2, "ip": 3, "component": 2},
            {"type": "stress", "element": 2, "ip": 3, "von_mises": True},
        ],
    }
    json.dump(request, open(os.path.join(_HERE, "sensitivity_request.json"), "w"), indent=2)
    shutil.copy(fx["base"], os.path.join(_HERE, "replay_record.json"))

    # ---- the actual Program-2 work: opaque binary + manifest + record + request
    result = run_request(fx["base"], fx["manifest"], request)
    json.dump(result, open(os.path.join(_HERE, "sensitivity_result.json"), "w"), indent=2)

    d = result["diagnostics"]
    print("=" * 74)
    print(" Program 2 (collaborator) replay sensitivity | model=%s drive=%s"
          % (result["model_id"], args.drive))
    print("=" * 74)
    print(" (JHU source never seen; only the opaque binary + manifest + record + request)")
    print(" contract version           : %s" % json.load(open(fx["manifest"]))["contract_version"])
    print(" replay stress vs production : max_rel = %.2e" % result["stress_check"]["max_rel"])
    print(" equilibrium ||R_free||      : %.2e" % result["equilibrium"]["free_norm"])
    print(" timings (s)                 : replay=%.4f solve=%.4f total=%.4f"
          % (d["timings_s"]["replay"], d["timings_s"]["solve"], d["timings_s"]["total"]))
    print("-" * 74)
    print(" dq/dp   (analytic  vs  central-FD of the collaborator's perturbed analyses)")
    print("-" * 74)
    recs = {k: json.load(open(fx[k])) for k in ("Ep", "Em", "nup", "num")}
    worst = 0.0
    for o in result["outputs"]:
        req = o["request"]
        tag = "vm" if req.get("von_mises") else req.get("component", "")
        for pn in ("E", "nu"):
            hp = fx["hE"] if pn == "E" else fx["hnu"]
            fd = (_q(recs["Ep" if pn == "E" else "nup"], req)
                  - _q(recs["Em" if pn == "E" else "num"], req)) / (2 * hp)
            an = o["dq_dp"][pn]
            rel = abs(an - fd) / max(abs(fd), 1e-30)
            worst = max(worst, rel if abs(fd) > 1e-3 else 0.0)
            print("   dq/d%-2s %-14s analytic=%+12.5e  FD=%+12.5e  %s"
                  % (pn, "%s[%s]" % (req["type"], tag), an, fd,
                     ("rel=%.1e" % rel) if abs(fd) > 1e-3 else "(~0 ok)"))
    print("-" * 74)
    print(" worst rel error on non-trivial sensitivities: %.2e -> %s"
          % (worst, "PASS" if worst < 1e-5 else "FAIL"))
    return 0 if worst < 1e-5 else 1


if __name__ == "__main__":
    sys.exit(main())
