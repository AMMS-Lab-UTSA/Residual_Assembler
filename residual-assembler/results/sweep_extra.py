#!/usr/bin/env python3
"""Compact 20-material summary table + a representative OTI-vs-FD plot for the
sweep. Reads results/figures/sweep_error_tables.json (produced by
build_sweep_tables.py) for the summary; marches one new material for the plot.

    UMAT_OTI_ROOT=~/Documents/UMAT_source_transformation python results/sweep_extra.py
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
import build_sweep_tables as B                       # noqa: E402
from residual_core.replay.path_material import PathMaterial, mises, dmises_dstress  # noqa: E402

FIG = os.path.join(RA, "results", "figures")
P1 = os.path.expanduser(os.environ.get("UMAT_OTI_ROOT", "~/Documents/UMAT_source_transformation"))


def _desc(mid):
    mdir = mid[5:] if mid.startswith("umat_") else mid
    return B.MODEL_DESC.get(mdir, mdir.replace("sweep_", " ").replace("_", " "))


def summary_table():
    data = json.load(open(os.path.join(FIG, "sweep_error_tables.json")))
    order = sorted(data, key=lambda m: (0 if m.startswith("umat_") else 1, m))
    rows = []
    for mid in order:
        d = data[mid]
        rows.append([mid, _desc(mid)[:38], str(len(d["program1"])),
                     "path-dep" if d.get("path_dependent") else "elastic",
                     "%.1e" % max(d["program1"].values()), "%.1e" % max(d["program2"].values())])

    def cc(v):
        v = float(v)
        return "#c9e7c0" if v < 1e-5 else ("#fde6b8" if v < 1e-3 else "#f4b8b8")
    colors = [["#ffffff"] * 4 + [cc(r[4]), cc(r[5])] for r in rows]
    fig, ax = plt.subplots(figsize=(13.2, 0.36 * len(rows) + 1.5)); ax.axis("off")
    hdr = ["Material", "Physics", "# params", "type", "P1 worst\n(DSIGMA_DP vs FD)", "P2 worst\n(residual vs FD)"]
    t = ax.table(cellText=rows, colLabels=hdr, cellColours=colors, colColours=[B.HDR] * 6,
                 loc="center", cellLoc="left")
    t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.5)
    for (rr, _c), cell in t.get_celld().items():
        if rr == 0:
            cell.set_text_props(color="white", fontweight="bold")
        cell.set_edgecolor("white")
    t.auto_set_column_width(list(range(6)))
    ng = sum(1 for r in rows if float(r[4]) < 1e-5 and float(r[5]) < 1e-5)
    ax.set_title("Whole framework validated across %d material models  (%d / %d exact to <1e-5 on both programs)"
                 % (len(rows), ng, len(rows)), fontsize=13.5, fontweight="bold", pad=14, color="#1a2a4a")
    ax.text(0.5, -0.03, "Program 1 = OTI DSIGMA_DP vs finite differences of the original UMAT.   "
            "Program 2 = residual-method d(sigma_vM)/dp vs full-analysis FD.   "
            "6 model materials + 14 authored/adapted (incl. real ICP ECO, ECL_TEMP, PCO).",
            transform=ax.transAxes, ha="center", fontsize=9, color="#5b6b7a", style="italic")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "table_sweep_summary.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    return len(rows), ng


def oti_vs_fd_plot(mdir="sweep_j2_kinematic", emax=0.02, N=100):
    MAT = os.path.join(P1, "oti_provider", "materials", mdir)
    c = json.load(open(os.path.join(MAT, "umat_%s_oti.json" % mdir)))
    obj = os.path.join(MAT, c["object"]["file"])
    props = json.load(open(os.path.join(MAT, "contract.json")))["validation"]["props_values"]
    names = [p["name"] for p in c["parameters"]]; pidx = [p["props_index"] - 1 for p in c["parameters"]]
    mat = PathMaterial(obj, c, workdir=tempfile.mkdtemp())
    dt = [1.0 / N] * N; de = [[emax / N, 0, 0, 0, 0, 0] for _ in range(N)]
    axial = np.cumsum([emax / N] * N) * 100
    oti = mat.march_oti(props, de, dt)
    W = np.array([props[i] for i in pidx])
    woti = np.zeros((N, len(names)))
    for k, o in enumerate(oti):
        g = dmises_dstress(o["stress"])
        for j in range(len(names)):
            woti[k, j] = (g @ o["dsigma_dp"][:, j]) * W[j]
    wfd = np.zeros((N, len(names)))
    for j, ip in enumerate(pidx):
        h = 1e-5 * (abs(props[ip]) or 1.0); pp = list(props); pp[ip] += h; pm = list(props); pm[ip] -= h
        rp = mat.march_regular(pp, de, dt); rm = mat.march_regular(pm, de, dt)
        wfd[:, j] = np.array([(mises(rp[n]["stress"]) - mises(rm[n]["stress"])) / (2 * h) for n in range(N)]) * W[j]
    idx = np.linspace(0, N - 1, 12).round().astype(int)
    cols = plt.cm.viridis(np.linspace(0.1, 0.85, len(names)))
    fig, axes = plt.subplots(1, len(names), figsize=(16, 4.2), sharex=True)
    worst = 0.0
    for j, nm in enumerate(names):
        a = axes[j]; a.plot(axial, woti[:, j], "-", color=cols[j], lw=3, label="OTI")
        a.plot(axial[idx], wfd[idx, j], "o", color="#222222", ms=7, mfc="white", mew=1.6, label="central FD")
        s = np.max(np.abs(wfd[:, j])); r = np.sqrt(np.mean((woti[:, j] - wfd[:, j]) ** 2)) / (s + 1e-30) if s > 0 else 0
        worst = max(worst, r)
        a.set_title(nm, fontsize=17, fontweight="bold"); a.tick_params(labelsize=13); a.grid(alpha=0.25)
        a.text(0.97, 0.05, "rel. RMSE %.1e" % r, transform=a.transAxes, ha="right", va="bottom",
               fontsize=13, fontweight="bold", color="#111827")
        a.set_xlabel("axial strain (%)", fontsize=14)
    axes[0].set_ylabel(r"$p\,\partial\sigma_{vM}/\partial p$  (MPa)", fontsize=15)
    axes[0].legend(loc="upper left", fontsize=13)
    fig.suptitle("OTI vs finite differences across the loading path  -  J2 plasticity with kinematic hardening "
                 "(new sweep material)   ·   worst rel. RMSE %.1e" % worst,
                 fontsize=16, fontweight="bold", color="#1a2a4a")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(FIG, "sweep_oti_vs_fd_j2kin.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    return worst


def main():
    n, ng = summary_table()
    w = oti_vs_fd_plot()
    print("summary table (%d materials, %d both-exact) + OTI-vs-FD plot (worst %.1e)" % (n, ng, w))
    return 0


if __name__ == "__main__":
    sys.exit(main())
