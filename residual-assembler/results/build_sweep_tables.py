#!/usr/bin/env python3
"""Expanded per-UMAT, per-parameter error tables for the ~20-UMAT sweep.

Auto-discovers every material under oti_provider/materials/ that has a completed
OTI contract (the 6 model materials + the sweep_* set + sweep_eco), and for each
runs BOTH validations exactly as build_error_tables.py does:

  Program 1 (transform):  OTI DSIGMA_DP vs central FD of the original UMAT.
  Program 2 (residual):   the residual method's d(sigma_vM)/dp vs FD of the full
                          controlled analysis (single C3D8, uniform deformation).

props0 comes from each transform contract's validation.props_values; path
dependence from the completed contract's history. A physics label is read from
results/figures/sweep_results.json when present, else the directory name.

    UMAT_OTI_ROOT=~/Documents/UMAT_source_transformation \
      python results/build_sweep_tables.py
"""
import glob
import json
import os
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, RA)
P1 = os.path.expanduser(os.environ.get("UMAT_OTI_ROOT", "~/Documents/UMAT_source_transformation"))
from residual_core.replay.path_material import PathMaterial, mises, dmises_dstress
from residual_core.replay.objlink import package_from_contract
from residual_core.replay.package import MaterialPackage
from residual_core.replay.abi import MaterialABI

FIGDIR = os.path.join(RA, "results", "figures"); os.makedirs(FIGDIR, exist_ok=True)
MATROOT = os.path.join(P1, "oti_provider", "materials")
HDR = "#4a6d8c"

# concise, accurate descriptions (the "real" ICP ones are self-contained adaptations,
# not the literal production source whose deep external call trees the transform can't take)
MODEL_DESC = {
    "m1_elastic": "Isotropic linear elasticity",
    "m2_cubic": "Cubic anisotropic elasticity",
    "m2_elastic3d": "Isotropic 3D elasticity",
    "m3_j2": "J2 plasticity, linear hardening",
    "m5_cpflow": "Thermally-activated viscoplastic flow",
    "m6_fcc": "FCC single crystal, 12 slip systems",
    "sweep_aniso_ortho": "Orthotropic elasticity",
    "sweep_lame_elastic": "Isotropic elasticity (Lame)",
    "sweep_transiso": "Transversely isotropic elasticity",
    "sweep_thermoelastic": "Thermo-elasticity",
    "sweep_damage_elastic": "Elasticity with scalar damage",
    "sweep_j2_kinematic": "J2 plasticity, kinematic hardening",
    "sweep_j2_combined": "J2 plasticity, combined hardening",
    "sweep_j2_bilinear": "J2 plasticity, bilinear hardening",
    "sweep_drucker_prager": "Drucker-Prager plasticity",
    "sweep_perzyna_linear": "Perzyna viscoplasticity",
    "sweep_maxwell_ve": "Maxwell viscoelasticity",
    "sweep_mooney_small": "Neo-Hookean (small strain)",
    "sweep_eco": "Cosserat elasticity (ICP ECO, adapted)",
    "sweep_real_ECL_TEMP": "Thermo-elasticity (ICP ECL_TEMP, adapted)",
    "sweep_real_PCO": "Couple-stress plasticity (ICP PCO, reformulated)",
}


def discover():
    """Return [(mdir, desc, props0, path_dep, completed_contract_path)] for every
    material that has both an OTI object and a completed contract."""
    meta = {}
    sr = os.path.join(FIGDIR, "sweep_results.json")
    if os.path.exists(sr):
        for r in json.load(open(sr)):
            meta[r.get("material_dir_name") or ("sweep_" + r["name"])] = r
    out = []
    for cpath in sorted(glob.glob(os.path.join(MATROOT, "*", "umat_*_oti.json"))):
        mdir = os.path.basename(os.path.dirname(cpath))
        try:
            c = json.load(open(cpath))
        except (OSError, ValueError):
            continue
        if c.get("schema") != "resasm_umat_oti_contract_v1":
            continue
        obj = os.path.join(os.path.dirname(cpath), c["object"]["file"])
        if not os.path.exists(obj):
            continue
        # props0: transform contract validation.props_values
        tc = os.path.join(MATROOT, mdir, "contract.json")
        props0 = None
        if os.path.exists(tc):
            props0 = json.load(open(tc)).get("validation", {}).get("props_values")
        if props0 is None and mdir in meta:
            props0 = meta[mdir].get("props_values")
        if props0 is None:
            continue
        path_dep = bool(c.get("history", {}).get("path_dependent")) or int(c["dimensions"].get("nstatev", 0)) > 0
        desc = MODEL_DESC.get(mdir) or (meta.get(mdir, {}).get("physics")) or mdir.replace("sweep_", "").replace("_", " ")
        out.append((mdir, desc, list(map(float, props0)), path_dep, cpath))
    return out


def _load(mdir, cpath):
    c = json.load(open(cpath)); obj = os.path.join(MATROOT, mdir, c["object"]["file"])
    names = [p["name"] for p in c["parameters"]]; pidx = [p["props_index"] - 1 for p in c["parameters"]]
    return c, obj, names, pidx


def prog1_path(mdir, cpath, props0):
    c, obj, names, pidx = _load(mdir, cpath)
    mat = PathMaterial(obj, c, workdir=tempfile.mkdtemp())
    N = 150; emax = 0.03; dt = [1.0 / N] * N
    dstran = [[emax / N, 0, 0, 0, 0, 0] for _ in range(N)]
    oti = mat.march_oti(props0, dstran, dt)
    dsig = {nm: np.array([o["dsigma_dp"][:, k] for o in oti]) for k, nm in enumerate(names)}
    fds, rmse = {}, {}
    for k, ip in enumerate(pidx):
        h = 1e-5 * (abs(props0[ip]) or 1.0)
        pp = list(props0); pp[ip] += h; pm = list(props0); pm[ip] -= h
        rp = mat.march_regular(pp, dstran, dt); rm = mat.march_regular(pm, dstran, dt)
        fd = np.array([(np.array(rp[n]["stress"]) - np.array(rm[n]["stress"])) / (2 * h) for n in range(N)])
        fds[names[k]] = fd; rmse[names[k]] = float(np.sqrt(np.mean((dsig[names[k]] - fd) ** 2)))
    gscale = max(np.max(np.abs(fds[nm])) for nm in names)
    return names, {nm: rmse[nm] / max(np.max(np.abs(fds[nm])), 1e-3 * gscale) for nm in names}


def prog1_stateless(mdir, cpath, props0):
    c, obj, names, pidx = _load(mdir, cpath)
    work = tempfile.mkdtemp(); base, manifest = package_from_contract(obj, c, work)
    pkg = MaterialPackage(manifest, base_dir=base); abi = MaterialABI(pkg.binary_path())
    seed = [int(p["props_index"]) for p in c["parameters"]]
    eps = np.array([1e-3, -3e-4, 2e-4, 1e-4, 5e-5, -2e-5])
    dsig = abi.eval_point(props0, seed, eps)["dstress_dseed"]
    mat = PathMaterial(obj, c, workdir=work)
    fds, aerr = {}, {}
    for k, ip in enumerate(pidx):
        h = 1e-6 * (abs(props0[ip]) or 1.0)
        pp = list(props0); pp[ip] += h; pm = list(props0); pm[ip] -= h
        sp = mat.march_regular(pp, [list(eps)], [1.0])[0]["stress"]
        sm = mat.march_regular(pm, [list(eps)], [1.0])[0]["stress"]
        fd = (np.array(sp) - np.array(sm)) / (2 * h)
        fds[names[k]] = fd; aerr[names[k]] = float(np.max(np.abs(dsig[:, k] - fd)))
    gscale = max(np.max(np.abs(fds[nm])) for nm in names)
    return names, {nm: aerr[nm] / max(np.max(np.abs(fds[nm])), 1e-3 * gscale) for nm in names}


def prog2_residual(mdir, cpath, props0, path_dep):
    c, obj, names, pidx = _load(mdir, cpath)
    if path_dep:
        mat = PathMaterial(obj, c, workdir=tempfile.mkdtemp())
        N = 120; emax = 0.03; dt = [1.0 / N] * N
        dstran = [[emax / N, 0, 0, 0, 0, 0] for _ in range(N)]
        oti = mat.march_oti(props0, dstran, dt)
        sig = oti[-1]["stress"]; dsig_dp = oti[-1]["dsigma_dp"]; g = dmises_dstress(sig)
        dsvm = {nm: float(g @ dsig_dp[:, k]) for k, nm in enumerate(names)}
        fds, aerr = {}, {}
        for k, ip in enumerate(pidx):
            h = 1e-5 * (abs(props0[ip]) or 1.0)
            pp = list(props0); pp[ip] += h; pm = list(props0); pm[ip] -= h
            rp = mat.march_regular(pp, dstran, dt); rm = mat.march_regular(pm, dstran, dt)
            fd = (mises(rp[-1]["stress"]) - mises(rm[-1]["stress"])) / (2 * h)
            fds[names[k]] = fd; aerr[names[k]] = abs(dsvm[names[k]] - fd)
        gscale = max(abs(fds[nm]) for nm in names) or 1.0
        return names, {nm: aerr[nm] / max(abs(fds[nm]), 1e-3 * gscale) for nm in names}
    work = tempfile.mkdtemp(); base, manifest = package_from_contract(obj, c, work)
    pkg = MaterialPackage(manifest, base_dir=base); abi = MaterialABI(pkg.binary_path())
    mat = PathMaterial(obj, c, workdir=work); seed = [int(p["props_index"]) for p in c["parameters"]]
    eps = np.array([2e-3, -6e-4, -6e-4, 0, 0, 0])
    res = abi.eval_point(props0, seed, eps); sig = res["stress"]; dsig = res["dstress_dseed"]; g = dmises_dstress(sig)
    fds, aerr = {}, {}
    for k, ip in enumerate(pidx):
        h = 1e-6 * (abs(props0[ip]) or 1.0)
        pp = list(props0); pp[ip] += h; pm = list(props0); pm[ip] -= h
        sp = mat.march_regular(pp, [list(eps)], [1.0])[0]["stress"]
        sm = mat.march_regular(pm, [list(eps)], [1.0])[0]["stress"]
        fd = (mises(sp) - mises(sm)) / (2 * h)
        fds[names[k]] = fd; aerr[names[k]] = abs(float(g @ dsig[:, k]) - fd)
    gscale = max(abs(fds[nm]) for nm in names) or 1.0
    return names, {nm: aerr[nm] / max(abs(fds[nm]), 1e-3 * gscale) for nm in names}


def _render(rows, title, fname, note):
    fig, ax = plt.subplots(figsize=(12.8, 0.30 * len(rows) + 1.5)); ax.axis("off")
    header = ["Material (file)", "Physics", "Parameter", "OTI vs FD (rel.)"]
    data, colors = [], []
    for mid, desc, param, err in rows:
        data.append([mid, desc[:34], param, "%.1e" % err if err > 0 else "0"])
        c = "#c9e7c0" if err < 1e-5 else ("#fde6b8" if err < 1e-3 else "#f4b8b8")
        colors.append(["#ffffff", "#ffffff", "#ffffff", c])
    tbl = ax.table(cellText=data, colLabels=header, cellColours=colors,
                   colColours=[HDR] * 4, loc="center", cellLoc="left")
    tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.28)
    for (r, cc), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_text_props(color="white", fontweight="bold")
        cell.set_edgecolor("white")
    tbl.auto_set_column_width([0, 1, 2, 3])
    ax.set_title(title, fontsize=12.5, fontweight="bold", pad=12, color="#1a2a4a")
    ax.text(0.5, -0.03, note, transform=ax.transAxes, ha="center", fontsize=9, color="#5b6b7a", style="italic")
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, fname), dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    mats = discover()
    print("discovered %d materials with OTI objects" % len(mats))
    p1_rows, p2_rows, summary = [], [], {}
    for mdir, desc, props0, path_dep, cpath in mats:
        mid = mdir if mdir.startswith("sweep_") else "umat_" + mdir
        try:
            n1, e1 = (prog1_path if path_dep else prog1_stateless)(mdir, cpath, props0)
            n2, e2 = prog2_residual(mdir, cpath, props0, path_dep)
        except Exception as ex:
            print("  %-22s SKIP (%s)" % (mdir, str(ex)[:60])); continue
        for nm in n1:
            p1_rows.append((mid, desc, nm, e1[nm])); p2_rows.append((mid, desc, nm, e2[nm]))
        summary[mid] = {"physics": desc, "path_dependent": path_dep, "program1": e1, "program2": e2}
        print("  %-22s P1 worst=%.1e  P2 worst=%.1e  (%d params)" % (mdir, max(e1.values()), max(e2.values()), len(n1)))
    json.dump(summary, open(os.path.join(FIGDIR, "sweep_error_tables.json"), "w"), indent=2)
    nmat = len(summary)
    w1 = max((max(v["program1"].values()) for v in summary.values()), default=0)
    w2 = max((max(v["program2"].values()) for v in summary.values()), default=0)
    _render(p1_rows, "Program 1  -  UMAT->OTI transformation across %d materials: DSIGMA_DP per parameter vs FD of the original UMAT  (worst %.1e)" % (nmat, w1),
            "table_program1_sweep.png",
            "Every differentiated parameter of every material, OTI stress derivative vs central finite differences of the separately compiled original UMAT.")
    _render(p2_rows, "Program 2  -  residual assembly across %d materials: d(sigma_vM)/dp per parameter vs FD of the full analysis  (worst %.1e)" % (nmat, w2),
            "table_program2_sweep.png",
            "Single C3D8, controlled deformation; the residual method assembles R, dR/dp, K, solves K du/dp = -dR/dp and evaluates d(sigma_vM)/dp vs full-analysis FD.")
    print("=> %d materials, P1 worst %.1e, P2 worst %.1e" % (nmat, w1, w2))
    print("tables -> table_program1_sweep.png, table_program2_sweep.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
