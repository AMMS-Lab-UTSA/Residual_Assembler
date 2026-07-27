#!/usr/bin/env python3
"""FCC single-crystal plasticity results (Program 1 + Program 2 material replay).

Marches the distributed OTI .obj for the 12-slip-system FCC crystal (m6_fcc) under
[100] uniaxial loading, extracting all TEN d(sigma_vM)/dp in one enriched OTI run
and validating each against central finite differences of the regular UMAT.
Produces, in the study's style:

  fcc_weighted_regimes.png   weighted sigma_vM sensitivities + sigma_vM, regime bands
  fcc_oti_vs_fd.png          9 non-trivial params, OTI vs FD, (p/s) ds/dp + rel RMSE
  fcc_timing.png             nominal / OTI / FD run cost + OTI breakdown
  fcc_results.json           per-parameter errors + timing (for the results tables)

    UMAT_OTI_ROOT=~/Documents/UMAT_source_transformation \
      python results/fcc_crystal_results.py
"""
import json, os, sys, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, RA)
P1 = os.path.expanduser(os.environ.get("UMAT_OTI_ROOT", "~/Documents/UMAT_source_transformation"))
from residual_core.replay.path_material import PathMaterial, mises, dmises_dstress

FIGDIR = os.path.join(RA, "results", "figures"); os.makedirs(FIGDIR, exist_ok=True)
MATDIR = os.path.join(P1, "oti_provider", "materials", "m6_fcc")

# [100] FCC crystal parameters (tuned to a reference-scale response).
PROPS = [168000.0, 121000.0, 75000.0, 13.0, 55.0, 800.0, 2.0, 1.4, 0.001, 0.05]
#         C11       C12       C44     g0    gsat  h0     a    q    gd0    m
G0_IDX, GSAT_IDX = 3, 4

# display order (like the reference) + LaTeX labels + slide-matched colours
DISPLAY = ["g0", "h0", "q", "gd0", "m", "gsat", "C11", "C12", "a"]   # C44 ~ 0, omitted
LATEX = {"g0": r"$g_0$", "h0": r"$h_0$", "q": r"$q$", "gd0": r"$\dot{\gamma}_0$",
         "m": r"$m$", "gsat": r"$g_{sat}$", "C11": r"$C_{11}$", "C12": r"$C_{12}$",
         "C44": r"$C_{44}$", "a": r"$a$"}
COLORS = {"g0": "#2a7fd4", "h0": "#5aa9e6", "q": "#f0a848", "gd0": "#f2d24b",
          "m": "#2fa66a", "gsat": "#7cc47c", "C11": "#d64b4b", "C12": "#e88b8b",
          "C44": "#999999", "a": "#8e5bd6"}


def main():
    contract = json.load(open(os.path.join(MATDIR, "umat_m6_fcc_oti.json")))
    obj = os.path.join(MATDIR, contract["object"]["file"])
    names = [p["name"] for p in contract["parameters"]]
    pidx = [p["props_index"] - 1 for p in contract["parameters"]]
    pval = np.array([PROPS[i] for i in pidx])
    mat = PathMaterial(obj, contract, workdir=os.path.join(FIGDIR, "_fcc_link"))
    nprm = mat.nparam

    N, emax = 250, 0.05
    dt = [1.0 / N] * N
    dstran = [[emax / N, 0, 0, 0, 0, 0] for _ in range(N)]
    strain = np.array([(n + 1) * emax / N for n in range(N)]) * 100.0   # %

    # ---- OTI: one EFFICIENT enriched march (state-alive, param-only OTI) --
    reg = mat.march_regular(PROPS, dstran, dt)                    # stress path (nominal work)
    fast = mat.march_fast(PROPS, dstran, dt)                      # all 10 sensitivities, one call
    reps = 10
    t0 = time.perf_counter()
    for _ in range(reps):
        mat.march_fast(PROPS, dstran, dt)
    t_oti = (time.perf_counter() - t0) / reps
    t0 = time.perf_counter()
    for _ in range(reps):
        mat.march_regular(PROPS, dstran, dt)
    t_reg = (time.perf_counter() - t0) / reps
    svm = np.array([mises(r["stress"]) for r in reg])
    dsvm_oti = np.array([dmises_dstress(reg[n]["stress"]) @ fast[n]["dsigma_dp"] for n in range(N)])
    gbar = np.array([r["statev"].mean() for r in reg])           # accumulated hardening
    dsvm_fd = np.zeros((N, nprm))
    for k, ip in enumerate(pidx):
        h = 1e-5 * (abs(PROPS[ip]) or 1.0)
        pp = list(PROPS); pp[ip] += h; pm = list(PROPS); pm[ip] -= h
        rp = mat.march_regular(pp, dstran, dt); rm = mat.march_regular(pm, dstran, dt)
        for n in range(N):
            dsvm_fd[n, k] = (mises(rp[n]["stress"]) - mises(rm[n]["stress"])) / (2 * h)

    idx = {nm: k for k, nm in enumerate(names)}
    rmse_rel = {}
    for nm in names:
        k = idx[nm]; sc = max(np.max(np.abs(dsvm_fd[:, k])), 1e-30)
        rmse_rel[nm] = float(np.sqrt(np.mean((dsvm_oti[:, k] - dsvm_fd[:, k]) ** 2)) / sc)

    # ---- regimes from hardening fraction (elastic/yield/hardening/saturation)
    hrange = PROPS[GSAT_IDX] - PROPS[G0_IDX]
    frac = gbar / hrange
    def first(th):
        w = np.where(frac > th)[0]; return strain[w[0]] if len(w) else strain[-1]
    b1, b2, b3 = first(0.005), first(0.05), first(0.30)     # yield, hardening, saturation
    regimes = [(strain[0], b1, "Elastic"), (b1, b2, "Yield"),
               (b2, b3, "Hardening"), (b3, strain[-1], "Saturation")]

    _plot_weighted(strain, svm, dsvm_oti, pval, names, idx, regimes)
    _plot_oti_vs_fd(strain, svm, dsvm_oti, dsvm_fd, pval, names, idx, rmse_rel)
    # fcc_timing.png is generated at the ANALYSIS level by cp_analysis_timing.py
    # (nominal analysis + enriched replay vs 2N+1 finite-difference analyses).

    summary = {
        "material": contract["model_id"], "orientation": "[100]", "n_slip_systems": 12,
        "parameters": names, "props": PROPS, "n_increments": N, "strain_max_pct": emax * 100,
        "sigma_vM_final_MPa": float(svm[-1]),
        "rmse_rel_oti_vs_fd": rmse_rel, "worst_rmse_rel": float(max(rmse_rel.values())),
        "regimes_pct": {lab: [float(a), float(b)] for a, b, lab in regimes},
        "cost": {"nominal_run_s": t_reg, "oti_run_s": t_oti,
                 "oti_vs_nominal": t_oti / max(t_reg, 1e-9), "fd_runs": 2 * nprm + 1,
                 "fd_total_vs_nominal": float(2 * nprm + 1)},
    }
    json.dump(summary, open(os.path.join(FIGDIR, "fcc_results.json"), "w"), indent=2)
    print("=" * 72)
    print(" FCC [100] single crystal — 10 sensitivities in one OTI run")
    print("=" * 72)
    print(" sigma_vM(5%%) = %.1f MPa ; regimes(%%): elastic<%.2f  yield<%.2f  hardening<%.2f  saturation"
          % (svm[-1], b1, b2, b3))
    for nm in DISPLAY + ["C44"]:
        print("  d(sigVM)/d%-4s OTI-vs-FD relRMSE = %.2e" % (nm, rmse_rel[nm]))
    print(" worst = %.2e" % max(rmse_rel.values()))
    print(" cost: nominal %.2f ms ; OTI %.2f ms (%.1fx) ; FD %d runs"
          % (t_reg * 1e3, t_oti * 1e3, t_oti / max(t_reg, 1e-9), 2 * nprm + 1))
    print(" figures + fcc_results.json ->", FIGDIR)
    return 0 if max(rmse_rel.values()) < 1e-3 else 1


def _plot_weighted(strain, svm, dsvm, pval, names, idx, regimes):
    weighted = np.zeros((len(strain), len(names)))
    for k in range(len(names)):
        weighted[:, k] = np.abs(dsvm[:, k]) * abs(pval[k])
    tot = np.maximum(weighted.sum(1, keepdims=True), 1e-30)
    wpct = 100.0 * weighted / tot
    order = [n for n in ["g0", "h0", "q", "gd0", "m", "gsat", "C11", "C12", "a", "C44"] if n in names]

    fig, axL = plt.subplots(figsize=(11.5, 4.8)); axR = axL.twinx()
    cum = np.zeros_like(strain)
    for nm in order:
        k = idx[nm]
        axL.fill_between(strain, cum, cum + wpct[:, k], color=COLORS[nm], alpha=0.85,
                         linewidth=0, label=LATEX[nm])
        cum = cum + wpct[:, k]
    axR.plot(strain, svm / 1000.0, color="black", lw=2.4, zorder=6)
    # regime bands + markers
    reg_bounds = [regimes[0][0]] + [r[1] for r in regimes]
    for i, (a, b, lab) in enumerate(regimes):
        if i > 0:
            axL.axvline(a, color="black", ls="--", lw=1.6, zorder=5)
        axL.text((a + b) / 2, 101.5, str(i + 1), ha="center", va="bottom",
                 fontsize=13, fontweight="bold", color="#333333")
    axL.set_ylim(0, 100); axR.set_ylim(0, max(svm) / 1000.0 * 1.05)
    axL.set_xlim(strain[0], strain[-1])
    axL.set_ylabel("Weighted Sensitivities (%)", fontweight="bold")
    axR.set_ylabel(r"$\sigma_{vM}$ (GPa)", fontweight="bold")
    axL.set_xlabel("Strain (%)", fontweight="bold")
    axL.set_title(r"FCC single crystal [100] : weighted von Mises $\sigma_{vM}$ sensitivities",
                  fontsize=13, fontweight="bold", color="#1a2a4a", pad=22)
    axL.legend(ncol=5, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.14), frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "fcc_weighted_regimes.png"), dpi=150,
                                    bbox_inches="tight"); plt.close(fig)


def _plot_oti_vs_fd(strain, svm, oti, fd, pval, names, idx, rmse_rel):
    fig, axes = plt.subplots(3, 3, figsize=(12, 8))
    for ax, nm in zip(axes.flat, DISPLAY):
        k = idx[nm]
        norm = pval[k] / np.maximum(svm, 1e-30)
        ax.plot(strain, norm * oti[:, k], color="#1f4e79", lw=2.4, label="OTI")
        ax.plot(strain, norm * fd[:, k], "--", color="#c02020", lw=1.8, label="FD")
        ax.set_title(LATEX[nm], fontsize=14)
        ax.set_xlabel("Strain (%)"); ax.set_ylabel(r"$(p/\sigma)\,\partial\sigma/\partial p$")
        ax.text(0.5, 0.08, "rel. RMSE = %.1e" % rmse_rel[nm], transform=ax.transAxes,
                ha="center", fontsize=11, fontweight="bold", color="#1a2a4a")
        ax.grid(alpha=0.3)
        if nm == "g0":
            ax.legend(fontsize=9, loc="upper right")
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "fcc_oti_vs_fd.png"), dpi=150)
    plt.close(fig)


def _plot_timing(t_reg, t_oti, nprm):
    nominal = t_reg * 1e3; oti = t_oti * 1e3
    n_fd = 2 * nprm + 1
    fd = n_fd * nominal
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw={"width_ratios": [1.4, 1]})
    # left: wall-time -- OTI (one enriched march) vs finite differences (2N+1 marches)
    bars = axL.bar(["Nominal march\n(no sensitivities)", "OTI residual method\n(all %d)" % nprm,
                    "Finite differences\n(all %d)" % nprm],
                   [nominal, oti, fd], color=["#8a8a8a", "#2a7fd4", "#c0202a"], width=0.62)
    axL.bar(["OTI residual method\n(all %d)" % nprm], [oti - nominal], bottom=[nominal],
            color="#2e8b57", width=0.62)
    for b, v in zip(bars, [nominal, oti, fd]):
        axL.text(b.get_x() + b.get_width() / 2, v + fd * 0.01, "%.0f ms" % v,
                 ha="center", va="bottom", fontsize=12, fontweight="bold")
    axL.set_ylabel("Time per run (ms)"); axL.set_ylim(0, fd * 1.15)
    axL.set_title("OTI gives all %d sensitivities in ONE march,\nabout %.1fx faster than finite differences"
                  % (nprm, fd / max(oti, 1e-9)), fontsize=12, fontweight="bold", color="#1a2a4a")
    # right: what the one OTI march spends time on
    axR.bar(["OTI march\n(%.0f ms)" % oti], [nominal], color="#8a8a8a", width=0.5, label="shared material work")
    axR.bar(["OTI march\n(%.0f ms)" % oti], [oti - nominal], bottom=[nominal], color="#2e8b57",
            width=0.5, label="extra for the %d sensitivities" % nprm)
    axR.text(0, nominal / 2, "%.0f ms" % nominal, ha="center", va="center", fontsize=11,
             fontweight="bold", color="white")
    axR.text(0, nominal + (oti - nominal) / 2, "%.0f ms" % (oti - nominal), ha="center",
             va="center", fontsize=11, fontweight="bold", color="white")
    axR.set_ylabel("Time per run (ms)")
    axR.set_title("What the OTI march spends time on", fontsize=12, fontweight="bold", color="#1a2a4a")
    axR.legend(fontsize=8.5, loc="upper left")
    ndir = 6 + int(nprm)
    axR.text(0.5, -0.16, "State-alive whole-path march: %d OTI directions (%d params + 6 strain),\n"
             "one call for all sensitivities + tangent. An optimized OTI library lowers this further."
             % (ndir, nprm), transform=axR.transAxes, ha="center", fontsize=8, color="#5b6b7a", style="italic")
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "fcc_timing.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
