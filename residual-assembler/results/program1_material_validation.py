#!/usr/bin/env python3
"""Program 1 (umat_oti) material-point validation figure: BOTH transformer
outputs vs finite differences, along the loading path.

The transformer returns, at every material point, the stress derivative
DSIGMA_DP and the state-variable derivative DSTATEV_DP. This figure validates
both against centered finite differences of the ORIGINAL UMAT, for the J2
linear-hardening model (M3): 4 parameters (E, nu, SIGY0, H) and one state
variable (equivalent plastic strain). OTI is drawn as lines, finite differences
as points -- exactly the OTI-vs-FD comparison style used in the deck.

The finite-difference reference is **Program 1's canonical one**
(`umat_oti.validation.fd_reference`): a parameter-scaled centered difference
whose step is chosen per parameter from the convergence of the FD sequence
itself, then Richardson-checked.  This script and the transformer's build-time
self-check therefore report the same quantity computed the same way -- they used
to each pick their own fixed step, which is how the FCC crystal model came to
"pass" here and "fail" there.

    python results/program1_material_validation.py
    UMAT_OTI_ROOT=/path/to/umat-oti python results/program1_material_validation.py
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, RA)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _p1_root import P1, fd_reference
fdref = fd_reference()
from residual_core.replay.path_material import PathMaterial, mises, dmises_dstress

FIGDIR = os.path.join(RA, "results", "figures"); os.makedirs(FIGDIR, exist_ok=True)
MATDIR = os.path.join(P1, "oti_provider", "materials", "m3_j2")
PVALS = {"E": 200000.0, "nu": 0.3, "SIGY0": 250.0, "H": 2000.0}
UNIT = {"E": "MPa", "nu": "-", "SIGY0": "MPa", "H": "MPa"}
N = 60
EPS_MAX = 0.010
FD_POINTS = 13


def make_props(contract, vals):
    props = [0.0] * int(contract["dimensions"]["nprops"])
    for p in contract["parameters"]:
        props[int(p["props_index"]) - 1] = vals[p["name"]]
    return props


def main():
    contract = json.load(open(os.path.join(MATDIR, "umat_m3_j2_oti.json")))
    obj = os.path.join(MATDIR, contract["object"]["file"])
    mat = PathMaterial(obj, contract, workdir=os.path.join(FIGDIR, "_m3_link"))
    names = [p["name"] for p in contract["parameters"]]
    props = make_props(contract, PVALS)

    # uniaxial-strain ramp (develops plasticity; J2 yield is pressure-independent)
    de = [0.0, 0.0, EPS_MAX / N, 0.0, 0.0, 0.0]
    path = [de] * N; dt = [1.0 / N] * N
    axial = np.cumsum([de[2]] * N) * 100.0                       # % axial strain

    # OTI march: dsigma_dp (6,nparam) and dstatev_dp (1,nparam) per increment
    oti = mat.march_oti(props, path, dt)
    svm = np.array([mises(o["stress"]) for o in oti])
    ep = np.array([o["statev"][0] for o in oti])
    dsvm_oti = np.zeros((N, len(names)))
    dep_oti = np.zeros((N, len(names)))
    for k, o in enumerate(oti):
        g = dmises_dstress(o["stress"])
        for j in range(len(names)):
            dsvm_oti[k, j] = g @ o["dsigma_dp"][:, j]
            dep_oti[k, j] = o["dstatev_dp"][0, j]

    # canonical converged centered finite differences of the ORIGINAL UMAT.
    # The step is selected per parameter from the FD sequence alone (never from
    # agreement with OTI) and Richardson-checked; see umat_oti.validation.fd_reference.
    def response(props_vec):
        reg = mat.march_regular(list(props_vec), path, dt)
        return {"svm": fdref.require_finite("original sigma_vM", [mises(r["stress"]) for r in reg]),
                "ep": fdref.require_finite("original EQPLAS", [r["statev"][0] for r in reg])}

    dsvm_fd = np.zeros((N, len(names))); dep_fd = np.zeros((N, len(names)))
    fd_evidence = []
    for j, nm in enumerate(names):
        idx_props = int(contract["parameters"][j]["props_index"])
        studies = fdref.parameter_studies(response, props, idx_props, parameter=nm)
        for block, target in (("svm", dsvm_fd), ("ep", dep_fd)):
            study = studies[block]
            fdref.apply_convergence_gate(study, 1e-6)
            if study.reference is None:
                raise RuntimeError("no converged FD reference for %s/%s: %s" % (nm, block, study.note))
            target[:, j] = study.reference
            fd_evidence.append(study.to_json())

    # parameter-weighted (dimensionless-ish) sensitivities: p * d(.)/dp
    W = np.array([PVALS[n] for n in names])
    wsvm_oti = dsvm_oti * W; wsvm_fd = dsvm_fd * W
    wep_oti = dep_oti * W;   wep_fd = dep_fd * W

    def relrmse(o, f):
        s = np.max(np.abs(f))
        return np.sqrt(np.mean((o - f) ** 2)) / (s + 1e-30) if s > 0 else 0.0

    idx = np.linspace(0, N - 1, FD_POINTS).round().astype(int)
    cols = plt.cm.viridis(np.linspace(0.1, 0.85, len(names)))
    fig, axes = plt.subplots(2, len(names), figsize=(18.0, 9.2), sharex=True)
    worst = 0.0

    def _rmse_label(ax, r):
        ax.text(0.98, 0.06, "rel. RMSE %.1e" % r, transform=ax.transAxes, fontsize=13,
                color="#111827", ha="right", va="bottom", fontweight="bold")

    for j, nm in enumerate(names):
        # top: DSIGMA_DP -> d(sigma_vM)/dp
        a = axes[0, j]
        a.plot(axial, wsvm_oti[:, j], "-", color=cols[j], lw=3.0, label="OTI")
        a.plot(axial[idx], wsvm_fd[idx, j], "o", color="#222222", ms=7, mfc="white", mew=1.6, label="central FD")
        r = relrmse(wsvm_oti[:, j], wsvm_fd[:, j]); worst = max(worst, r)
        a.set_title("%s   (%s)" % (nm, UNIT[nm]), fontsize=18, fontweight="bold")
        _rmse_label(a, r)
        a.grid(alpha=0.25); a.tick_params(labelsize=14)
        # bottom: DSTATEV_DP -> d(eq. plastic strain)/dp
        b = axes[1, j]
        b.plot(axial, wep_oti[:, j], "-", color=cols[j], lw=3.0)
        b.plot(axial[idx], wep_fd[idx, j], "o", color="#222222", ms=7, mfc="white", mew=1.6)
        r2 = relrmse(wep_oti[:, j], wep_fd[:, j]); worst = max(worst, r2)
        _rmse_label(b, r2)
        b.set_xlabel("axial strain  (%)", fontsize=16)
        b.grid(alpha=0.25); b.tick_params(labelsize=14)
    axes[0, 0].set_ylabel("DSIGMA_DP  (stress)\n$p\\,\\partial\\sigma_{vM}/\\partial p$   (MPa)", fontsize=16)
    axes[1, 0].set_ylabel("DSTATEV_DP  (state)\n$p\\,\\partial(\\mathrm{EQPLAS})/\\partial p$", fontsize=16)
    axes[0, 0].legend(loc="upper left", fontsize=14, framealpha=0.9)
    fig.suptitle("Program 1 (umat_oti) — both transformer outputs vs finite differences    ·    worst rel. RMSE %.1e"
                 % worst, fontsize=18, fontweight="bold", color="#1a2a4a", y=0.985)
    fig.text(0.5, 0.945,
             "J2 plasticity (M3):   normal variable = von Mises stress $\\sigma_{vM}$ (DSIGMA_DP)      ·      "
             "state variable = accumulated equivalent plastic strain, EQPLAS (DSTATEV_DP)      ·      "
             "OTI lines vs central-FD points",
             ha="center", va="top", fontsize=13.5, color="#333333")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = os.path.join(FIGDIR, "program1_dsigma_dstatev_vs_fd.png")
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
    import importlib
    prov = importlib.import_module("umat_oti.runtime.provenance")
    metric_def = {
        "worst_rel_rmse": {
            "quantity": "OTI-vs-FD relative RMSE of parameter-weighted DSIGMA_DP and DSTATEV_DP",
            "formula": "RMSE_over_path(p*dOTI/dp - p*dFD/dp) / max(|p*dFD/dp|)  (per parameter, worst reported)",
            "aggregation": "root-mean-square over the loading path; worst over parameters and blocks",
            "path": "%d-increment uniaxial-strain ramp to %.1f%%" % (N, EPS_MAX * 100),
            "normalization": "max |parameter-weighted FD| over the path",
            "terminology": "OTI-vs-FD relative discrepancy (RMSE), not error vs an analytical derivative",
        },
        "fd_reference": "umat_oti.validation.fd_reference (converged per-parameter step, Richardson-checked)",
    }
    json.dump({"material": "m3_j2", "params": names, "worst_rel_rmse": float(worst),
               "sigy0": PVALS["SIGY0"], "svm_final": float(svm[-1]), "ep_final": float(ep[-1]),
               "fd_method": "umat_oti.validation.fd_reference (converged, per-parameter step)",
               "method_version": "adaptive-fd-v1",
               "metric_definition": metric_def,
               "provenance": {
                   "git_commits": {"umat_oti": prov.git_commit(P1),
                                   "residual_assembler": prov.git_commit(RA)},
                   "environment": prov.build_environment(),
                   "props_values": PVALS,
                   "fd_ladder": list(fdref.DEFAULT_LADDER),
                   "timestamp_utc": prov.timestamp_utc()},
               "fd_evidence": fd_evidence},
              open(os.path.join(FIGDIR, "program1_material_validation.json"), "w"), indent=2)
    print("worst rel. RMSE (DSIGMA_DP & DSTATEV_DP vs FD) = %.2e" % worst)
    print("sigma_vM final = %.1f MPa,  eq. plastic strain final = %.4f" % (svm[-1], ep[-1]))
    print("finite-difference reference (umat_oti.validation.fd_reference):")
    print("  %-8s %-5s %13s %8s %10s %12s" % ("param", "block", "h tested", "sel. h",
                                              "FD self", "Richardson"))
    for e in fd_evidence:
        hs = e["h_tested"]
        print("  %-8s %-5s %13s %8.0e %10.1e %12.1e"
              % (e["parameter"], e["block"], "%.0e..%.0e" % (hs[0], hs[-1]),
                 e["selected_h_rel"], e["convergence_delta"], e["richardson_gap"]))
    print("figure ->", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
