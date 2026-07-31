#!/usr/bin/env python3
"""Exercise the Program 2 sensitivity_job.json CLI end-to-end and FD-check it.

  1. Program 1 (subprocess) builds umat_m1_elastic_oti.obj + completed json.
  2. Build a converged record (regular UMAT) with node/element sets, write to disk.
  3. Write a real resasm_sensitivity_job_v1 job with region+reduction outputs.
  4. run_job -> reduced DQ_DP for nodal-displacement(avg), reaction(sum),
     stress(volume_average).
  5. FD-check every reduced sensitivity against perturbed regular-UMAT records.

    UMAT_OTI_ROOT=~/Documents/UMAT_source_transformation python validation/run_sensitivity_job.py
"""
import json, os, subprocess, sys, tempfile
import numpy as np

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, RA); sys.path.insert(0, os.path.dirname(__file__))
P1 = os.path.expanduser(os.environ.get("UMAT_OTI_ROOT", "~/Documents/UMAT_source_transformation"))

import reference_analysis as ref
from residual_core.replay.job import run_job

MAT, PN, DRIVE = "m1_elastic", ["E", "nu"], "disp"


def _reduce_fd(recs, param, h, kind, members_reduction):
    """Central-difference the SAME reduction the job applies, from regular-UMAT records."""
    def per_member(rec):
        vals = []
        for req in members_reduction["reqs"]:
            vals.append(ref.response_from_record(rec, req))
        return np.array(vals)
    vp, vm = per_member(recs[param + "+"]), per_member(recs[param + "-"])
    dv = (vp - vm) / (2 * h)
    w = members_reduction.get("weights")
    if members_reduction["reduction"] == "sum":
        return float(np.sum(dv))
    if members_reduction["reduction"] == "volume_average":
        return float(np.dot(w, dv) / np.sum(w))
    return float(np.mean(dv))


def main():
    cdir = os.path.join(P1, "oti_provider", "materials", MAT)
    contract_in = json.load(open(os.path.join(cdir, "contract.json")))
    props0 = list(contract_in["validation"]["props_values"])
    r = subprocess.run([sys.executable, os.path.join(P1, "oti_provider", "umat_transform.py"),
                        "build", os.path.join(cdir, "contract.json"),
                        "--props", ",".join(str(x) for x in props0)],
                       capture_output=True, text=True)
    obj = os.path.join(cdir, contract_in["output"]["object"])
    cjson = os.path.join(cdir, contract_in["output"]["contract"])
    if r.returncode != 0 or not os.path.exists(obj):
        print("Program 1 build FAILED:\n", (r.stdout + r.stderr)[-800:]); return 1
    completed = json.load(open(cjson))
    d = completed["dimensions"]

    work = tempfile.mkdtemp(prefix="job_demo_")
    reg = ref.RegularUmat(_link_regular(obj, completed, work), d["ntens"], d["nprops"], d["nstatev"])
    mid = completed["model_id"]; rhash = completed.get("regular_source_hash")
    m, dm, f, nodes, elems = ref.cube_model(DRIVE)
    steps = {p["name"]: 1e-6 * (abs(props0[p["props_index"]-1]) or 1.0) for p in completed["parameters"]}
    recs = {"base": ref.make_record(m, dm, reg, props0, f, nodes, elems, DRIVE, mid, rhash)}
    for p in completed["parameters"]:
        j = p["props_index"]-1; h = steps[p["name"]]
        pp = list(props0); pp[j] += h; pm = list(props0); pm[j] -= h
        recs[p["name"]+"+"] = ref.make_record(m, dm, reg, pp, f, nodes, elems, DRIVE, mid, rhash)
        recs[p["name"]+"-"] = ref.make_record(m, dm, reg, pm, f, nodes, elems, DRIVE, mid, rhash)
    rec_path = os.path.join(work, "cube.resrec.json"); json.dump(recs["base"], open(rec_path, "w"))

    job = {
        "schema": "resasm_sensitivity_job_v1",
        "material": {"oti_umat": obj, "contract": cjson},
        "analysis": {"record": rec_path},
        "location": {"step": "Step-1", "increment": "last"},
        "parameters": [{"name": p} for p in PN],
        "outputs": [
            {"output_id": "top_uz_avg", "type": "nodal_displacement",
             "region": {"node_set": "LOADED_FACE"}, "component": "U3", "reduction": "average"},
            {"output_id": "base_reaction_sum", "type": "reaction_force",
             "region": {"node_set": "FIXED_FACE"}, "component": "RF3", "reduction": "sum"},
            {"output_id": "domain_s33_vavg", "type": "integration_point_stress",
             "region": {"element_set": "DOMAIN"}, "component": "S33",
             "integration_points": "all", "reduction": "volume_average"}],
        "result": {"summary_json": True, "summary_csv": True, "output_dir": os.path.join(work, "out")},
    }
    job_path = os.path.join(work, "sensitivity_job.json"); json.dump(job, open(job_path, "w"), indent=2)

    result = run_job(job_path)

    # FD plans mirroring the three reductions
    plans = {
        "top_uz_avg": {"reduction": "average",
                       "reqs": [{"type": "displacement", "node": n, "dof": 3} for n in (9, 10, 11, 12)]},
        "base_reaction_sum": {"reduction": "sum",
                              "reqs": [{"type": "reaction", "node": n, "dof": 3} for n in (1, 2, 3, 4)]},
    }
    # stress volume-average FD plan needs detJ*w weights
    from residual_core.formulations import c3d8_kernel as kern
    model = ref.cube_model(DRIVE)[0]
    reqs, wts = [], []
    for eid in (1, 2):
        Xe = model.coords_of(model.elements[eid].connectivity)
        for k in range(8):
            reqs.append({"type": "stress", "element": eid, "ip": k, "component": 2})
            wts.append(kern.b_matrix_reference(Xe, kern.ABAQUS_C3D8_GAUSS.points[k])[1]
                       * kern.ABAQUS_C3D8_GAUSS.weights[k])
    plans["domain_s33_vavg"] = {"reduction": "volume_average", "reqs": reqs, "weights": np.array(wts)}

    print("=" * 78); print(" PROGRAM 2  sensitivity_job.json  (reduced DQ_DP vs regular-UMAT FD)")
    print("=" * 78)
    worst = 0.0
    for o in result["outputs"]:
        plan = plans[o["output_id"]]
        line = "  %-20s value=%+.5e" % (o["output_id"], o["value"])
        for p in PN:
            h = steps[p]
            fd = _reduce_fd(recs, p, h, o["output_id"], plan)
            got = o["dq_dp"][p]
            rel = abs(got - fd) / max(1.0, abs(fd))
            worst = max(worst, rel)
            line += "  d/d%s=%+.4e(fd %+.4e rel %.1e)" % (p, got, fd, rel)
        print(line)
    print("-" * 78)
    print(" ||R_free||=%.2e  replay_stress_max_rel=%.2e  worst_rel=%.2e -> %s"
          % (result["equilibrium"]["free_norm"], result["stress_check"]["max_rel"],
             worst, "PASS" if worst < 1e-4 else "FAIL"))
    print(" result written:", os.path.join(work, "out", "sensitivity_result.json"))
    return 0 if worst < 1e-4 else 1


def _link_regular(obj, completed, work):
    from residual_core.replay.objlink import package_from_contract
    base, manifest = package_from_contract(obj, completed, work)
    return os.path.join(base, manifest["binaries"]["oti"]["path"])


if __name__ == "__main__":
    sys.exit(main())
