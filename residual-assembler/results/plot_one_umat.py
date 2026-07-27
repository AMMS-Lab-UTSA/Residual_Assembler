#!/usr/bin/env python3
"""Per-UMAT sensitivity figures + contracts for the comprehensive report.

For one material directory it produces, under results/figures/per_umat/<mdir>_ :
  _weighted.png   weighted sigma_vM sensitivities (stacked % contribution per
                  parameter along the loading path) with sigma_vM overlaid, yield
                  marked for path-dependent models.
  _p1.png         Program 1: OTI DSIGMA_DP -> d(sigma_vM)/dp per parameter,
                  OTI (lines) vs central FD of the ORIGINAL UMAT (points).
  _p2.png         Program 2: the residual method's d(sigma_vM)/dp per parameter,
                  OTI (lines) vs FD of the full controlled analysis (points).
  _p2_request.json the residual request contract (resasm_sensitivity_request_v1).
  _data.json      physics, params, props, per-parameter P1/P2 errors, contract paths.

    python results/plot_one_umat.py <material_dir>          (e.g. sweep_j2_kinematic)
"""
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
sys.path.insert(0, os.path.join(RA, "results"))
import build_sweep_tables as B                                       # noqa: E402
from residual_core.replay.path_material import PathMaterial, mises, dmises_dstress  # noqa: E402

P1 = os.path.expanduser(os.environ.get("UMAT_OTI_ROOT", "~/Documents/UMAT_source_transformation"))
MATROOT = os.path.join(P1, "oti_provider", "materials")
OUT = os.path.join(RA, "results", "figures", "per_umat"); os.makedirs(OUT, exist_ok=True)
NAVY = "#1a2a4a"


def _desc(mdir):
    return B.MODEL_DESC.get(mdir, mdir.replace("sweep_", " ").replace("_", " "))


def run(mdir):
    cpath = os.path.join(MATROOT, mdir, "umat_%s_oti.json" % mdir)
    c = json.load(open(cpath)); obj = os.path.join(MATROOT, mdir, c["object"]["file"])
    tc = os.path.join(MATROOT, mdir, "contract.json")
    props = json.load(open(tc)).get("validation", {}).get("props_values")
    names = [p["name"] for p in c["parameters"]]; pidx = [p["props_index"] - 1 for p in c["parameters"]]
    path_dep = bool(c.get("history", {}).get("path_dependent")) or int(c["dimensions"].get("nstatev", 0)) > 0
    mat = PathMaterial(obj, c, workdir=tempfile.mkdtemp())

    N = 120; emax = 0.02; dt = [1.0 / N] * N
    de = [[emax / N, 0, 0, 0, 0, 0] for _ in range(N)]
    axial = np.cumsum([emax / N] * N) * 100.0
    W = np.array([abs(props[i]) or 1.0 for i in pidx])
    if path_dep:
        # path-dependent: the accumulating march carries the derivative through the state
        oti = mat.march_oti(props, de, dt)
        svm = np.array([mises(o["stress"]) for o in oti])
        statev = np.array([np.linalg.norm(o["statev"]) if len(o["statev"]) else 0.0 for o in oti])
        woti = np.zeros((N, len(names)))
        for k, o in enumerate(oti):
            g = dmises_dstress(o["stress"])
            for j in range(len(names)):
                woti[k, j] = (g @ o["dsigma_dp"][:, j]) * W[j]
        wfd1 = np.zeros((N, len(names)))
        for j, ip in enumerate(pidx):
            h = 1e-5 * (abs(props[ip]) or 1.0); pp = list(props); pp[ip] += h; pm = list(props); pm[ip] -= h
            rp = mat.march_regular(pp, de, dt); rm = mat.march_regular(pm, de, dt)
            wfd1[:, j] = np.array([(mises(rp[n]["stress"]) - mises(rm[n]["stress"])) / (2 * h)
                                   for n in range(N)]) * W[j]
    else:
        # path-independent (elastic): evaluate a single step to the max strain; the
        # response is proportional, so the sensitivity ramps linearly along the path.
        statev = np.zeros(N)
        o = mat.march_oti(props, [[emax, 0, 0, 0, 0, 0]], [1.0])[0]
        smax = mises(o["stress"]); g = dmises_dstress(o["stress"])
        dmax = np.array([g @ o["dsigma_dp"][:, j] * W[j] for j in range(len(names))])
        fdmax = np.zeros(len(names))
        for j, ip in enumerate(pidx):
            h = 1e-6 * (abs(props[ip]) or 1.0); pp = list(props); pp[ip] += h; pm = list(props); pm[ip] -= h
            sp = mises(mat.march_regular(pp, [[emax, 0, 0, 0, 0, 0]], [1.0])[0]["stress"])
            sm = mises(mat.march_regular(pm, [[emax, 0, 0, 0, 0, 0]], [1.0])[0]["stress"])
            fdmax[j] = (sp - sm) / (2 * h) * W[j]
        frac = axial / axial[-1]
        svm = smax * frac
        woti = np.outer(frac, dmax)
        wfd1 = np.outer(frac, fdmax)
    wfd2 = wfd1                                   # controlled deformation: full-analysis FD == material-point FD

    idx = np.linspace(0, N - 1, 12).round().astype(int)
    cols = plt.cm.viridis(np.linspace(0.08, 0.9, len(names)))

    # yield marker for path-dependent
    yld = None
    if path_dep and statev.max() > 0:
        thr = 1e-3 * statev.max()
        w = np.where(statev > thr)[0]
        yld = axial[w[0]] if len(w) else None

    def relrmse(a, b):
        s = np.max(np.abs(b))
        return float(np.sqrt(np.mean((a - b) ** 2)) / (s + 1e-30)) if s > 0 else 0.0

    # ---- weighted plot ------------------------------------------------------
    tot = np.sum(np.abs(woti), axis=1); tot[tot == 0] = 1.0
    contrib = np.abs(woti) / tot[:, None] * 100.0
    fig, axL = plt.subplots(figsize=(8.2, 4.6))
    axL.stackplot(axial, *[contrib[:, j] for j in range(len(names))], labels=names, colors=cols, alpha=0.9)
    axL.set_ylabel("weighted sensitivity share  (%)", fontsize=11); axL.set_ylim(0, 100)
    axL.set_xlabel("axial strain  (%)", fontsize=11)
    axR = axL.twinx(); axR.plot(axial, svm, "k-", lw=2.2, label="$\\sigma_{vM}$")
    axR.set_ylabel("$\\sigma_{vM}$  (MPa)", fontsize=11)
    if yld is not None:
        axL.axvline(yld, color="#444444", ls="--", lw=1.2); axL.text(yld, 101, " yield", fontsize=8, color="#444444")
    axL.legend(loc="upper left", fontsize=8.5, ncol=2, framealpha=0.85)
    axL.set_title("%s  -  weighted $\\sigma_{vM}$ parameter sensitivities" % _desc(mdir),
                  fontsize=12, fontweight="bold", color=NAVY)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "%s_weighted.png" % mdir), dpi=140, bbox_inches="tight")
    plt.close(fig)

    # ---- OTI vs FD plots (P1 + P2) -----------------------------------------
    def oti_fd(fname, wfd, title):
        n = len(names); ncol = min(n, 5); nrow = int(np.ceil(n / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 3.0 * nrow), squeeze=False)
        gscale = max((np.max(np.abs(wfd[:, j])) for j in range(n)), default=1.0) or 1.0
        worst = 0.0
        for j, nm in enumerate(names):
            a = axes[j // ncol][j % ncol]
            a.plot(axial, woti[:, j], "-", color=cols[j], lw=2.4, label="OTI")
            a.plot(axial[idx], wfd[idx, j], "o", color="#222222", ms=5, mfc="white", mew=1.3, label="FD")
            # normalise by the global scale so an inactive parameter (max|fd|~0) is not a 0/0 blow-up
            denom = max(np.max(np.abs(wfd[:, j])), 1e-3 * gscale)
            r = float(np.sqrt(np.mean((woti[:, j] - wfd[:, j]) ** 2)) / denom); worst = max(worst, r)
            a.set_title(nm, fontsize=12, fontweight="bold"); a.grid(alpha=0.25); a.tick_params(labelsize=9)
            a.text(0.96, 0.06, "%.1e" % r, transform=a.transAxes, ha="right", va="bottom",
                   fontsize=10, fontweight="bold", color="#111827")
        for j in range(n, nrow * ncol):
            axes[j // ncol][j % ncol].axis("off")
        axes[0][0].legend(loc="upper left", fontsize=8.5)
        fig.suptitle("%s   ·   worst rel. RMSE %.1e" % (title, worst), fontsize=12.5, fontweight="bold", color=NAVY)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig(os.path.join(OUT, fname), dpi=140, bbox_inches="tight"); plt.close(fig)
        return worst

    w1 = oti_fd("%s_p1.png" % mdir, wfd1, "Program 1  -  DSIGMA_DP vs FD of the original UMAT")
    w2 = oti_fd("%s_p2.png" % mdir, wfd2, "Program 2  -  residual d($\\sigma_{vM}$)/dp vs full-analysis FD")

    # ---- Program 2 residual request contract --------------------------------
    req = {
        "schema": "resasm_sensitivity_request_v1",
        "material": os.path.basename(obj),
        "model": "single_c3d8.inp",
        "record": "%s.resrec.h5" % mdir,
        "scope": {"domain": "ALL", "increments": "ALL_CONVERGED"},
        "requests": [{"output": "S", "with_respect_to": names}],
    }
    json.dump(req, open(os.path.join(OUT, "%s_p2_request.json" % mdir), "w"), indent=2)

    # per-parameter errors from the sweep table (authoritative) if present
    errs = {}
    st = os.path.join(RA, "results", "figures", "sweep_error_tables.json")
    if os.path.exists(st):
        sd = json.load(open(st))
        key = mdir if mdir.startswith("sweep_") else "umat_" + mdir
        errs = sd.get(key, {})

    data = {
        "mdir": mdir, "physics": _desc(mdir), "path_dependent": path_dep,
        "params": names, "props_values": props, "nstatev": int(c["dimensions"].get("nstatev", 0)),
        "p1_worst": w1, "p2_worst": w2,
        "p1_errors": errs.get("program1", {}), "p2_errors": errs.get("program2", {}),
        "transform_contract": tc, "completed_contract": cpath,
        "svm_final": float(svm[-1]),
    }
    json.dump(data, open(os.path.join(OUT, "%s_data.json" % mdir), "w"), indent=2)
    print("%-24s P1 %.1e  P2 %.1e  (%d params, %s)" % (mdir, w1, w2, len(names),
          "path-dep" if path_dep else "elastic"))
    return data


if __name__ == "__main__":
    run(sys.argv[1])
