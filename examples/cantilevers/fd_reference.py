#!/usr/bin/env python3
"""Central-difference reference from the perturbed Abaqus reruns.

    python fd_reference.py j2 [WORK] -> WORK/j2/fd_reference.json, .csv
                                        (WORK: where run_fd.sh ran; default .)

For every parameter p and relative step h with both reruns completed, the
central difference  dq/dp ~ (q(p(1+h)) - q(p(1-h))) / (2 h p)  is taken for:

  tip_RF2[step]        total reaction force on the tip face (N)
  U2_mid_top[step]     vertical displacement at the free node x=L/2, top, z=0
  S11_root_ip[step]    S11 at element 1 (root, bottom), integration point 1
  vm_root_max[step]    max von Mises over the root element column
                       (a max is not differentiable where the arg-max switches;
                        reported, but flagged as such)

A value is "resolved" when the three step sizes agree (plateau): the spread of
the h=0.01 and h=0.005 estimates relative to their magnitude is recorded, and
the h=0.02 estimate shows whether truncation error is still visible. Field
output in an ODB is single precision; the floor that imposes is estimated from
the float32 unit roundoff times the output magnitude divided by (2 h p).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gen_cantilever import MODELS  # noqa: E402


def outputs(npz):
    d = np.load(npz)
    coords = d["coords"]
    L = coords[:, 0].max(); H = coords[:, 1].max()
    tip = np.isclose(coords[:, 0], L)
    mid_top = np.where(np.isclose(coords[:, 0], L / 2) & np.isclose(coords[:, 1], H)
                       & np.isclose(coords[:, 2], 0))[0][0]
    s = d["S"]
    vm = np.sqrt(0.5 * ((s[..., 0] - s[..., 1]) ** 2 + (s[..., 1] - s[..., 2]) ** 2
                        + (s[..., 2] - s[..., 0]) ** 2) + 3 * (s[..., 3] ** 2 + s[..., 4] ** 2 + s[..., 5] ** 2))
    conn = d["conn"]; xs = coords[conn - 1][:, :, 0].mean(axis=1)
    root_col = np.isclose(xs, xs.min())
    return {
        "tip_RF2": d["RF"][:, tip, 1].sum(axis=1),
        "U2_mid_top": d["U"][:, mid_top, 1],
        "S11_root_ip": s[:, 0, 0, 0],
        "vm_root_max": vm[:, root_col, :].reshape(vm.shape[0], -1).max(axis=1),
    }


def main(model, work="."):
    here = Path(work).resolve()
    spec = MODELS[model]
    base = outputs(here / model / f"cantilever_{model}_nominal_fields.npz")
    rows = list(csv.DictReader(open(here / model / "fd_manifest.tsv"), delimiter="\t"))
    done = {(r["parameter"], r["sign"], r["rel_step"]): r for r in rows if r["status"] == "completed"}
    eps32 = np.finfo(np.float32).eps
    result = {"model": model, "props": spec["props"], "outputs": {}, "note": __doc__.split("\n\n")[2]}
    table = []
    for pi, name in enumerate(spec["prop_names"]):
        p0 = spec["props"][pi]
        per_h = {}
        for h in ("0.02", "0.01", "0.005"):
            if (name, "p", h) in done and (name, "m", h) in done:
                qp = outputs(here / model / f"{done[(name, 'p', h)]['job']}_fields.npz")
                qm = outputs(here / model / f"{done[(name, 'm', h)]['job']}_fields.npz")
                per_h[h] = {k: (qp[k] - qm[k]) / (2 * float(h) * p0) for k in qp}
        if not per_h:
            continue
        for key in base:
            series = {h: v[key].tolist() for h, v in per_h.items()}
            floor = (eps32 * np.abs(base[key]) / (2 * 0.005 * p0)).tolist()
            result["outputs"].setdefault(key, {})[name] = {"by_step_size": series, "float32_floor": floor}
            if "0.01" in per_h and "0.005" in per_h:
                a, b = per_h["0.01"][key], per_h["0.005"][key]
                for step in range(1, len(a)):
                    mag = max(abs(a[step]), abs(b[step]), 1e-300)
                    table.append({"output": key, "parameter": name, "step": step,
                                  "fd_h0.01": a[step], "fd_h0.005": b[step],
                                  "fd_h0.02": per_h.get("0.02", {}).get(key, [np.nan] * len(a))[step],
                                  "plateau_rel_spread": abs(a[step] - b[step]) / mag,
                                  "float32_floor": floor[step]})
    (here / model / "fd_reference.json").write_text(json.dumps(result, indent=1))
    with open(here / model / "fd_reference.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0]) if table else ["output"])
        w.writeheader(); w.writerows(table)
    print(f"{model}: {len(table)} FD rows written; parameters: {sorted({r['parameter'] for r in table})}")


if __name__ == "__main__":
    main(sys.argv[1], *sys.argv[2:3])
