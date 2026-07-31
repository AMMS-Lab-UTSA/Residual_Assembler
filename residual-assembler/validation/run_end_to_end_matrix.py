#!/usr/bin/env python3
"""End-to-end two-program validation matrix (corrected .obj architecture).

Only the distributed .obj + completed contract cross the boundary:
  Program 1 (subprocess): build umat_<name>_oti.obj + umat_<name>_oti.json.
  Program 2 (import):     LINK the .obj (objlink), build reference base/perturbed
                          analyses by calling the REGULAR UMAT (reg_eval), REPLAY
                          via UMAT_OTI_EVAL to get DU_DP/DQ_DP, and finite-
                          difference-validate against the regular-UMAT analyses.

No isotropic_D or hard-coded constitutive law anywhere; no ZIP; no Program-1
import. Abaqus is not required for these offline reference analyses.

    UMAT_OTI_ROOT=~/Documents/UMAT_source_transformation python validation/run_end_to_end_matrix.py
"""
import json, os, subprocess, sys, tempfile
import numpy as np

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, RA); sys.path.insert(0, os.path.dirname(__file__))
P1 = os.path.expanduser(os.environ.get("UMAT_OTI_ROOT", "~/Documents/UMAT_source_transformation"))

from residual_core.replay import replay_sensitivities, MaterialPackage, ReplayRecord
from residual_core.replay.objlink import package_from_contract
import reference_analysis as ref

MATERIALS = [
    ("V1_M1_elastic_C3D8", "m1_elastic", ["E", "nu"]),
    ("V2_M2_cubic_C3D8",   "m2_cubic",   ["C11", "C12", "C44"]),
]
OUTPUTS = [{"type": "displacement", "node": 11, "dof": 1},
           {"type": "reaction", "node": 11, "dof": 3},
           {"type": "stress", "element": 2, "ip": 3, "component": 2}]


def run_one(vid, matname, pnames, drive):
    cdir = os.path.join(P1, "oti_provider", "materials", matname)
    contract = json.load(open(os.path.join(cdir, "contract.json")))
    props0 = list(contract["validation"]["props_values"])
    # ---- Program 1: build the .obj (subprocess) ----
    r = subprocess.run([sys.executable, os.path.join(P1, "oti_provider", "umat_transform.py"),
                        "build", os.path.join(cdir, "contract.json"),
                        "--props", ",".join(str(x) for x in props0)],
                       capture_output=True, text=True)
    obj = os.path.join(cdir, contract["output"]["object"])
    cjson = os.path.join(cdir, contract["output"]["contract"])
    if r.returncode != 0 or not os.path.exists(obj):
        return {"pairing": vid, "final_status": "BLOCKED",
                "note": "Program 1 build failed", "log": (r.stdout + r.stderr)[-500:]}
    completed = json.load(open(cjson))
    # ---- Program 2: LINK the .obj + build reference analyses via regular UMAT ----
    tmp = tempfile.mkdtemp(prefix="e2e_")
    base_dir, manifest = package_from_contract(obj, completed, tmp)
    d = completed["dimensions"]
    reg = ref.RegularUmat(os.path.join(base_dir, manifest["binaries"]["oti"]["path"]),
                          d["ntens"], d["nprops"], d["nstatev"])
    mid, rhash = manifest["model_id"], manifest["binaries"]["regular"]["hash"]
    m, dm, f, nodes, elems = ref.cube_model(drive)
    steps = {p["name"]: 1e-6 * (abs(props0[p["props_index"]-1]) or 1.0)
             for p in completed["parameters"]}
    recs = {"base": ref.make_record(m, dm, reg, props0, f, nodes, elems, drive, mid, rhash)}
    for p in completed["parameters"]:
        j = p["props_index"] - 1; h = steps[p["name"]]
        pp = list(props0); pp[j] += h; pm = list(props0); pm[j] -= h
        recs[p["name"] + "+"] = ref.make_record(m, dm, reg, pp, f, nodes, elems, drive, mid, rhash)
        recs[p["name"] + "-"] = ref.make_record(m, dm, reg, pm, f, nodes, elems, drive, mid, rhash)
    # ---- Program 2: replay -> DU_DP / DQ_DP ----
    res = replay_sensitivities(ReplayRecord(recs["base"]), MaterialPackage(manifest, base_dir=base_dir),
                               pnames, OUTPUTS)
    worst = 0.0
    for p in completed["parameters"]:
        h = steps[p["name"]]
        du_fd = (np.array(recs[p["name"]+"+"]["increments"][0]["u"])
                 - np.array(recs[p["name"]+"-"]["increments"][0]["u"])) / (2*h)
        worst = max(worst, np.max(np.abs(res.du(p["name"]) - du_fd)) / max(1.0, np.max(np.abs(du_fd))))
    for o in res.outputs:
        for p in completed["parameters"]:
            h = steps[p["name"]]
            fd = (ref.response_from_record(recs[p["name"]+"+"], o["request"])
                  - ref.response_from_record(recs[p["name"]+"-"], o["request"])) / (2*h)
            if abs(fd) > 1e-3:
                worst = max(worst, abs(o["dq_dp"][p["name"]] - fd) / abs(fd))
    return {"pairing": vid, "material": mid, "element": "C3D8_small", "drive": drive,
            "path_dependence": bool(d["nstatev"] > 0),
            "replay_stress_max_rel": res.stress_check["max_rel"],
            "equilibrium_free_norm": res.equilibrium_free_norm,
            "worst_error": float(worst), "final_status": "PASS" if worst < 1e-4 else "FAIL"}


if __name__ == "__main__":
    rows = []
    for vid, mat, pn in MATERIALS:
        for drive in ("disp", "force"):
            rows.append(run_one("%s_%s" % (vid, drive), mat, pn, drive))
    print("=" * 92)
    print(" END-TO-END MATRIX  (.obj + completed contract only; regular UMAT for FD; no Abaqus)")
    print("=" * 92)
    for r in rows:
        if r["final_status"] == "BLOCKED":
            print(" %-26s BLOCKED: %s" % (r["pairing"], r.get("note")))
        else:
            print(" %-26s %-16s %-5s replayS=%.1e ||Rf||=%.1e worst=%.2e -> %s" % (
                r["pairing"], r["material"], r["drive"], r["replay_stress_max_rel"],
                r["equilibrium_free_norm"], r["worst_error"], r["final_status"]))
    os.makedirs(os.path.join(RA, "validation", "out"), exist_ok=True)
    json.dump(rows, open(os.path.join(RA, "validation", "out", "matrix.json"), "w"), indent=2)
    sys.exit(0 if all(r["final_status"] == "PASS" for r in rows) else 1)
