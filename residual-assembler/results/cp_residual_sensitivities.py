#!/usr/bin/env python3
"""Residual-method parameter sensitivities for the thermally-activated crystal-
plasticity flow model (M5), single integration point under simple shear.

Consumes ONLY the distributed OTI .obj + completed contract (Program-1 output).
Marches the transformed UMAT over a shear path E12: 0 -> 0.05 and, in ONE enriched
run, extracts d(sigma_vM)/dp for the six flow parameters (tau0, dG, q, p, gam0, H)
via OTI.  Central finite differences of the REGULAR UMAT over the same path are
the independent reference.  Produces three figures in the study's style:

  1. weighted_sigvm_sensitivities.png  -- stacked weighted %-sensitivities + sigma_vM
  2. dsigvm_dp_oti_vs_fd.png           -- per-parameter OTI vs FD with RMSE/max
  3. normalized_cost.png               -- one OTI run vs 2*6+1 FD runs

    UMAT_OTI_ROOT=~/Documents/UMAT_source_transformation \
      python results/cp_residual_sensitivities.py
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
MATDIR = os.path.join(P1, "oti_provider", "materials", "m5_cpflow")

# legend order + slide-matched colours: tau0, dG, q, p, gam0, H
LABELS = [r"$\tau_0$", r"$\Delta G$", r"$q$", r"$p$", r"$\dot{\gamma}_0$", r"$H$"]
COLORS = ["#8e5bd6", "#2a7fd4", "#d64b4b", "#f2c230", "#f08a3c", "#2fa6a0"]
DESCR = ["Friction stress", "Activation energy", "Glide exponent q",
         "Glide exponent p", "Reference slip rate", "Hardening modulus"]


def build_path(n=50, e12_max=0.05):
    """Simple-shear increments (engineering gamma12 = 2*E12), constant dt."""
    dgam = (2.0 * e12_max) / n
    dstran = [[0.0, 0.0, 0.0, dgam, 0.0, 0.0] for _ in range(n)]
    dt = [1.0] * n
    e12 = np.array([(k + 1) * dgam / 2.0 for k in range(n)])
    return dstran, dt, e12


def main():
    contract = json.load(open(os.path.join(MATDIR, "umat_m5_cpflow_oti.json")))
    obj = os.path.join(MATDIR, contract["object"]["file"])
    props0 = list(contract["validation"]["props_values"])
    pidx = [p["props_index"] - 1 for p in contract["parameters"]]
    pval = np.array([props0[i] for i in pidx])
    mat = PathMaterial(obj, contract)
    dstran, dt, e12 = build_path()
    nprm = mat.nparam

    # ---- OTI: one EFFICIENT state-alive march -> sigma_vM and d sigma_vM/dp -
    reg = mat.march_regular(props0, dstran, dt)    # stress path (nominal work)
    fast = mat.march_fast(props0, dstran, dt)      # all sensitivities, one call
    oti = reg                                      # (stress source for the curve)
    reps = 20
    t0 = time.perf_counter()
    for _ in range(reps):
        mat.march_fast(props0, dstran, dt)
    t_oti = (time.perf_counter() - t0) / reps
    t0 = time.perf_counter()
    for _ in range(reps):
        mat.march_regular(props0, dstran, dt)
    t_reg = (time.perf_counter() - t0) / reps
    svm = np.array([mises(s["stress"]) for s in reg])
    dsvm_oti = np.zeros((len(reg), nprm))
    for n in range(len(reg)):
        g = dmises_dstress(reg[n]["stress"])
        dsvm_oti[n] = g @ fast[n]["dsigma_dp"]     # chain rule: dσvM/dp

    # ---- FD reference: 2*nprm regular marches -----------------------------
    dsvm_fd = np.zeros((len(oti), nprm))
    for k, ip in enumerate(pidx):
        h = 1e-4 * (abs(props0[ip]) or 1.0)
        pp = list(props0); pp[ip] += h; pm = list(props0); pm[ip] -= h
        rp = mat.march_regular(pp, dstran, dt); rm = mat.march_regular(pm, dstran, dt)
        for n in range(len(oti)):
            dsvm_fd[n, k] = (mises(rp[n]["stress"]) - mises(rm[n]["stress"])) / (2 * h)

    rmse = np.sqrt(np.mean((dsvm_oti - dsvm_fd) ** 2, axis=0))
    scale = np.maximum(np.max(np.abs(dsvm_fd), axis=0), 1e-30)
    rmse_rel = rmse / scale

    # weighted (dimensionless) sensitivities: |dσvM/dp * p|
    weighted = np.abs(dsvm_oti) * np.abs(pval)[None, :]
    tot = np.maximum(weighted.sum(axis=1, keepdims=True), 1e-30)
    wpct = 100.0 * weighted / tot

    _plot_weighted(e12, svm, wpct)
    _plot_oti_vs_fd(e12, dsvm_oti, dsvm_fd, rmse_rel)
    _plot_cost(t_oti, t_reg, nprm)

    summary = {
        "material": contract["model_id"], "n_increments": len(oti),
        "e12_max": float(e12[-1]), "parameters": mat.params,
        "sigma_vM_MPa_final": float(svm[-1]),
        "rmse_rel_oti_vs_fd": {mat.params[k]: float(rmse_rel[k]) for k in range(nprm)},
        "worst_rmse_rel": float(rmse_rel.max()),
        "cost": {"oti_run_s": t_oti, "regular_run_s": t_reg,
                 "oti_vs_regular": t_oti / max(t_reg, 1e-9),
                 "fd_runs": 2 * nprm + 1},
    }
    json.dump(summary, open(os.path.join(FIGDIR, "cp_sensitivity_summary.json"), "w"), indent=2)
    print("=" * 70)
    print(" CP residual-method sensitivities (single IP, simple shear)")
    print("=" * 70)
    print(" sigma_vM(E12=%.3f) = %.1f MPa" % (e12[-1], svm[-1]))
    for k in range(nprm):
        print("  d(sigVM)/d%-5s  OTI vs FD rel RMSE = %.2e" % (mat.params[k], rmse_rel[k]))
    print(" worst OTI-vs-FD rel RMSE = %.2e" % rmse_rel.max())
    print(" cost: OTI 1 enriched run = %.3fs ; regular run = %.3fs ; ratio %.2fx ; FD needs %d runs"
          % (t_oti, t_reg, t_oti / max(t_reg, 1e-9), 2 * nprm + 1))
    print(" figures ->", FIGDIR)
    return 0 if rmse_rel.max() < 1e-3 else 1


def _plot_weighted(e12, svm, wpct):
    fig, axL = plt.subplots(figsize=(7.2, 4.6))
    axR = axL.twinx()
    # stacked weighted-% areas on the right axis
    cum = np.zeros_like(e12)
    for k in range(wpct.shape[1]):
        axR.fill_between(e12, cum, cum + wpct[:, k], color=COLORS[k], alpha=0.85,
                         linewidth=0, label=LABELS[k])
        cum = cum + wpct[:, k]
    axL.plot(e12, svm / 1000.0, color="#7a1f1f", lw=2.4, zorder=5)
    axL.set_xlabel(r"Shear strain $E_{12}$")
    axL.set_ylabel(r"Stress $\sigma_{vM}$ (GPa)")
    axR.set_ylabel("Weighted sensitivities (%)")
    axR.set_ylim(0, 100); axL.set_xlim(e12[0], e12[-1])
    axL.set_title(r"Weighted $\sigma_{vM}$ Parameter Sensitivities (Single Element)")
    axR.legend(ncol=3, fontsize=8, loc="lower center", framealpha=0.9,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "weighted_sigvm_sensitivities.png"), dpi=150)
    plt.close(fig)


def _plot_oti_vs_fd(e12, oti, fd, rmse_rel):
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.2))
    for k, ax in enumerate(axes.flat):
        ax.plot(e12, oti[:, k], color="#1f5fd4", lw=2.0, label="OTI")
        ax.plot(e12, fd[:, k], "o", color="#d62728", ms=3.2, mfc="none", label="FD")
        ax.set_title(LABELS[k], fontsize=13)
        ax.set_xlabel(r"$E_{12}$"); ax.set_ylabel(r"$d\sigma_{vM}/dp$")
        ax.text(0.05, 0.86, "RMSE/max: %.1e" % rmse_rel[k], transform=ax.transAxes,
                fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
        ax.legend(fontsize=8, loc="lower left")
        ax.grid(alpha=0.25)
    fig.suptitle("Crystal-plasticity flow: OTI vs central finite differences", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(os.path.join(FIGDIR, "dsigvm_dp_oti_vs_fd.png"), dpi=150); plt.close(fig)


def _plot_cost(t_oti, t_reg, nprm):
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    oti_x = t_oti / max(t_reg, 1e-9)
    fd_x = float(2 * nprm + 1)                      # 13 single runs
    bars = ax.bar(["Nominal\n(no sens.)", "OTI\n(1 march, all %d)" % nprm,
                   "FD Total\n(%d runs)" % (2 * nprm + 1)],
                  [1.0, oti_x, fd_x], color=["#8a8a8a", "#2a7fd4", "#d64b4b"], width=0.62)
    ax.set_ylabel("Normalized cost (relative to a plain UMAT run)")
    ax.set_title("OTI: all %d sensitivities in one march, %.1f× faster than FD" % (nprm, fd_x / max(oti_x, 1e-9)))
    for b, v in zip(bars, [1.0, oti_x, fd_x]):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.2, u"×%.1f" % v,
                ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax.set_ylim(0, fd_x * 1.18)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "normalized_cost.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
