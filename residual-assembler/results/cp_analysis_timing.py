#!/usr/bin/env python3
"""Residual-method timing at the ANALYSIS level (the honest, reference-style view).

A nominal FE analysis solves the nonlinear system incrementally with Newton
iterations. Finite-difference sensitivities re-run that whole analysis 2*N+1
times. The OTI residual method instead REPLAYS the converged record once (an
enriched material pass, no re-solve) and does a single linear sensitivity solve.

So the cost decomposes as, per the study's figure:
  nominal run          : one nonlinear analysis (Newton)               (shared)
  OTI residual method  : shared work + the enriched-pass overhead      (~a few x)
  finite differences   : (2N+1) nominal analyses

This times a single C3D8 crystal-plasticity element under axial force control and
renders the two-panel normalized-timing figure.

    UMAT_OTI_ROOT=~/Documents/UMAT_source_transformation \
      python results/cp_analysis_timing.py
"""
import json, os, sys, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, RA)
P1 = os.path.expanduser(os.environ.get("UMAT_OTI_ROOT", "~/Documents/UMAT_source_transformation"))
from residual_core.formulations import c3d8_kernel as kern
from residual_core.replay.path_material import PathMaterial, mises, dmises_dstress

FIGDIR = os.path.join(RA, "results", "figures"); os.makedirs(FIGDIR, exist_ok=True)
MATDIR = os.path.join(P1, "oti_provider", "materials", "m6_fcc")
PROPS = [168000.0, 121000.0, 75000.0, 13.0, 55.0, 800.0, 2.0, 1.4, 0.001, 0.05]

XE = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
               [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
FIXED = [0, 1, 2, 3]          # base nodes z=0 -> encastre
TOP = [4, 5, 6, 7]            # top nodes z=1 -> axial force


def _dofs():
    fixed = []
    for n in FIXED:
        fixed += [3 * n, 3 * n + 1, 3 * n + 2]
    free = [d for d in range(24) if d not in fixed]
    return np.array(free, int)


def nominal_analysis(mat, props, N, Fmax, dt):
    """Incremental Newton solve (regular UMAT). Returns the per-IP strain path,
    the converged u path, and the average Newton iteration count."""
    pts, wts = kern.ABAQUS_C3D8_GAUSS.points, kern.ABAQUS_C3D8_GAUSS.weights
    B = [kern.b_matrix_reference(XE, pts[k]) for k in range(8)]     # (B, detJ)
    free = _dofs()
    u = np.zeros(24); nt = mat.ntens; ns = mat.nstatev
    sig = np.zeros((8, nt)); sv = np.zeros((8, ns))
    fext_dir = np.zeros(24)
    for n in TOP:
        fext_dir[3 * n + 2] = 0.25                                  # axial force split
    eps_path = np.zeros((8, N, nt)); newton_iters = []
    for inc in range(N):
        u_prev = u.copy(); sig_prev = sig.copy(); sv_prev = sv.copy()
        Fext = fext_dir * (Fmax * (inc + 1) / N)
        for it in range(30):
            Fint = np.zeros(24); K = np.zeros((24, 24))
            for k in range(8):
                Bk, detJ = B[k]
                dstran = Bk @ (u - u_prev)
                s, svk, dd = mat._reg_step(props, sig_prev[k], sv_prev[k], dstran, dt)
                sig[k] = s; sv[k] = svk
                Fint += Bk.T @ s * (detJ * wts[k])
                K += Bk.T @ dd @ Bk * (detJ * wts[k])
            R = Fint - Fext
            if np.linalg.norm(R[free]) < 1e-6 * (np.linalg.norm(Fext) + 1e-9):
                break
            du = np.linalg.solve(K[np.ix_(free, free)], -R[free]); u[free] += du
        newton_iters.append(it + 1)
        for k in range(8):
            eps_path[k, inc] = kern.b_matrix_reference(XE, pts[k])[0] @ u
    return eps_path, np.mean(newton_iters)


def main():
    contract = json.load(open(os.path.join(MATDIR, "umat_m6_fcc_oti.json")))
    obj = os.path.join(MATDIR, contract["object"]["file"])
    mat = PathMaterial(obj, contract, workdir=os.path.join(FIGDIR, "_fcc_link"))
    nprm = mat.nparam
    N = 40; dt = 1.0 / N
    pts, wts = kern.ABAQUS_C3D8_GAUSS.points, kern.ABAQUS_C3D8_GAUSS.weights
    B = [kern.b_matrix_reference(XE, pts[k]) for k in range(8)]

    # ---- nominal analysis (measure time + Newton iterations) --------------
    reps = 5
    t0 = time.perf_counter()
    for _ in range(reps):
        eps_path, nnewt = nominal_analysis(mat, PROPS, N, 200.0, dt)
    t_nominal = (time.perf_counter() - t0) / reps
    # per-IP incremental strain path (for replay)
    dstran_ip = []
    for k in range(8):
        ep = eps_path[k]
        d = np.vstack([ep[0], np.diff(ep, axis=0)])
        dstran_ip.append([list(row) for row in d])
    dtl = [dt] * N

    # ---- OTI residual method: replay (enriched) + assemble + sensitivity solve
    t0 = time.perf_counter()
    for _ in range(reps):
        # replay every IP once (enriched, all sensitivities + tangent)
        fields = [mat.march_fast(PROPS, dstran_ip[k], dtl) for k in range(8)]
        # assemble K and dR/dp at the converged (last) increment
        free = _dofs()
        K = np.zeros((24, 24)); dR = np.zeros((24, nprm))
        for k in range(8):
            Bk, detJ = B[k]
            K += Bk.T @ fields[k][-1]["ddsdde"] @ Bk * (detJ * wts[k])
            dR[:, :] += (Bk.T @ fields[k][-1]["dsigma_dp"]) * (detJ * wts[k])
        dudp = np.zeros((24, nprm))
        dudp[free, :] = np.linalg.solve(K[np.ix_(free, free)], -dR[free, :])
    t_oti = (time.perf_counter() - t0) / reps

    # OTI residual method = the primal analysis you already have (shared) + the
    # enriched replay + linear sensitivity solve (extra). FD = (2N+1) analyses.
    t_extra = t_oti
    t_oti_total = t_nominal + t_extra
    n_fd = 2 * nprm + 1
    t_fd = n_fd * t_nominal
    _plot(t_nominal, t_extra, t_fd, nprm, nnewt)

    summary = {"n_increments": N, "avg_newton_iters": float(nnewt), "n_params": nprm,
               "nominal_ms": t_nominal * 1e3, "extra_ms": t_extra * 1e3,
               "oti_total_ms": t_oti_total * 1e3, "fd_total_ms": t_fd * 1e3,
               "oti_vs_nominal": t_oti_total / t_nominal, "fd_vs_oti": t_fd / t_oti_total,
               "fd_runs": n_fd}
    json.dump(summary, open(os.path.join(FIGDIR, "cp_analysis_timing.json"), "w"), indent=2)
    print("=" * 68)
    print(" CP single-element analysis-level timing (%d increments, avg %.1f Newton its)" % (N, nnewt))
    print("=" * 68)
    print(" nominal analysis          = %.1f ms" % (t_nominal * 1e3))
    print(" OTI residual method       = %.1f ms  (shared %.1f + extra %.1f; %.1fx nominal)"
          % (t_oti_total * 1e3, t_nominal * 1e3, t_extra * 1e3, t_oti_total / t_nominal))
    print(" finite differences        = %.1f ms  (%d analyses)" % (t_fd * 1e3, n_fd))
    print(" OTI is %.1fx faster than finite differences" % (t_fd / t_oti_total))
    print(" figure -> fcc_timing.png")
    return 0


def _plot(t_nom, t_extra, t_fd, nprm, nnewt):
    # normalize every time to the nominal run (implementation independent)
    nominal = 1.0
    extra = t_extra / t_nom
    shared = 1.0
    oti = shared + extra
    fd = t_fd / t_nom
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.6, 5.0), gridspec_kw={"width_ratios": [1.45, 1]})
    # left: nominal / OTI / FD (OTI stacked: shared + extra), normalized
    labels = ["Nominal run\n(no sensitivities)", "OTI residual method\n(all %d)" % nprm,
              "Finite differences\n(all %d)" % nprm]
    axL.bar(labels, [nominal, shared, fd], color=["#8a8a8a", "#8a8a8a", "#c0202a"], width=0.62)
    axL.bar([labels[1]], [extra], bottom=[shared], color="#2e8b57", width=0.62)
    for x, v in zip(range(3), [nominal, oti, fd]):
        axL.text(x, v + fd * 0.012, u"×%.1f" % v, ha="center", va="bottom", fontsize=13, fontweight="bold")
    axL.set_ylabel("Execution time (× nominal run)"); axL.set_ylim(0, fd * 1.15)
    axL.set_title("OTI gives all %d sensitivities in one run,\nabout %.0fx faster than finite differences"
                  % (nprm, fd / max(oti, 1e-9)), fontsize=12.5, fontweight="bold", color="#1a2a4a")
    # right: what the OTI run spends time on (normalized)
    axR.bar([u"OTI run\n(×%.1f)" % oti], [shared], color="#8a8a8a", width=0.5)
    axR.bar([u"OTI run\n(×%.1f)" % oti], [extra], bottom=[shared], color="#2e8b57", width=0.5)
    axR.text(0, shared / 2, u"×%.1f" % shared, ha="center", va="center", fontsize=12, fontweight="bold", color="white")
    axR.text(0, shared + extra / 2, u"×%.2f" % extra, ha="center", va="center", fontsize=12, fontweight="bold", color="white")
    axR.annotate("extra for the\n%d sensitivities\n(×%.2f)" % (nprm, extra),
                 xy=(0.28, shared + extra / 2), xytext=(0.62, shared + extra * 0.6),
                 fontsize=11, fontweight="bold", color="#1a5e2a", ha="left", va="center")
    axR.annotate("shared\n(any run\ndoes this)", xy=(0.28, shared / 2), xytext=(0.62, shared * 0.5),
                 fontsize=10, color="#5b6b7a", ha="left", va="center")
    axR.set_ylabel("Execution time (× nominal run)"); axR.set_xlim(-0.6, 1.1)
    axR.set_ylim(0, oti * 1.18)
    axR.set_title("What the OTI run spends time on", fontsize=12.5, fontweight="bold", color="#1a2a4a")
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "fcc_timing.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
