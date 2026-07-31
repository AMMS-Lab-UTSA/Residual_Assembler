#!/usr/bin/env python3
"""Expand + validate the supported element formulations (deck next-step #2).

The residual framework was C3D8-only. `solid3d_kernel` adds the common Abaqus 3D
solid family -- reduced integration (C3D8R, C3D20R), quadratic hexes (C3D20),
and tetrahedra (C3D4, C3D10) -- all reusing the SAME 3D OTI UMAT and the SAME
residual assembly. This script VALIDATES each element type two ways:

  A. Homogeneous patch test (all 6 elements). A uniform strain is prescribed on
     every node, so d(u)/dp = 0 and the assembled reaction-force sensitivity must
     equal the material-point derivative integrated over the element. Compares
     OTI d(reaction)/dp against centered finite differences of the material.
     This validates the element ASSEMBLY (B matrix, integration rule, Jacobian,
     Voigt bookkeeping) x the OTI material derivative for every element type,
     including the reduced-integration ones (no solve -> no hourglass issue).

  B. Full residual solve (fully-integrated elements). A single element is loaded
     under force control so there are free DOFs; the method assembles K and
     dR/dp, solves K du/dp = -dR/dp, and evaluates d(tip displacement)/dp. This
     exercises the K^-1 solve and is compared against full-analysis finite
     differences. Reduced-integration single elements are rank-deficient (a
     standard hourglass property), so they are validated at the assembly level
     (test A) and require hourglass control / a mesh for a standalone solve.

    UMAT_OTI_ROOT=~/Documents/UMAT_source_transformation \
      python results/element_formulations.py
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, RA)
P1 = os.path.expanduser(os.environ.get("UMAT_OTI_ROOT", "~/Documents/UMAT_source_transformation"))
from residual_core.formulations import solid3d_kernel as k3
from residual_core.replay.path_material import PathMaterial

FIGDIR = os.path.join(RA, "results", "figures"); os.makedirs(FIGDIR, exist_ok=True)
ELEMS = ["C3D8", "C3D8R", "C3D20", "C3D20R", "C3D4", "C3D10"]
FULLY_INTEGRATED = {"C3D8", "C3D20", "C3D4", "C3D10"}

M2_PROPS = {"C11": 168000.0, "C12": 121000.0, "C44": 75000.0}   # copper-like cubic elastic
# a homogeneous deformation with shear so C11, C12 AND C44 are all active
FGRAD = np.array([[0.0020, 0.0008, -0.0012],
                  [0.0005, 0.0015, 0.0007],
                  [0.0010, -0.0006, 0.0025]])
LOAD_DIR = np.array([0.25, 0.15, 1.0]); LOAD_DIR = LOAD_DIR / np.linalg.norm(LOAD_DIR)
F_TOTAL = 90.0


def ref_coords(etype, distort=0.06):
    e = k3.element(etype)
    if e.nnode in (8, 20):
        nat = k3._HEX_CORNERS if e.nnode == 8 else k3.HEX20_NODES
        Xe = 0.5 * (nat + 1.0)
    else:
        Xe = (k3.TET4_NODES if e.nnode == 4 else k3.TET10_NODES).copy()
    rng = np.random.RandomState(3)
    d = distort * rng.randn(*Xe.shape)
    d[Xe[:, 2] <= Xe[:, 2].min() + 1e-9] = 0.0                    # keep zmin face planar
    return Xe + d


def make_props(contract, values):
    props = [0.0] * int(contract["dimensions"]["nprops"])
    for p in contract["parameters"]:
        props[int(p["props_index"]) - 1] = values[p["name"]]
    return props


def replay_fields(mat, props, dstran_path, dt_path):
    if getattr(mat, "has_march", False):
        return mat.march_fast(props, dstran_path, dt_path)
    return mat.march_oti(props, dstran_path, dt_path)


# ===========================================================================
# Test A: homogeneous patch test -- assembled reaction sensitivity vs FD
# ===========================================================================
def patch_test(mat, contract, etype):
    Xe = ref_coords(etype)
    e = k3.element(etype); nn = e.nnode; ndof = 3 * nn
    pts, wts = k3.gauss_rule(etype); nip = len(wts)
    Bs = [k3.b_matrix(Xe, pts[k], etype) for k in range(nip)]
    u = (Xe @ FGRAD.T).reshape(-1)                               # homogeneous field
    eps_ip = [Bs[k][0] @ u for k in range(nip)]
    # reaction weight: bottom-face (z=min) dofs projected onto LOAD_DIR
    z = Xe[:, 2]; bottom = np.where(z <= z.min() + 1e-6)[0]
    w = np.zeros(ndof)
    for n in bottom:
        w[3 * n:3 * n + 3] = LOAD_DIR
    names = [p["name"] for p in contract["parameters"]]
    props = make_props(contract, M2_PROPS)
    # OTI: d(reaction)/dp = (dFint/dp)^T w, dFint/dp = sum_k B^T dsigma/dp detJ w_k
    dFdp = np.zeros((ndof, len(names)))
    for k in range(nip):
        f = replay_fields(mat, props, [list(eps_ip[k])], [1.0])[-1]
        Bk, detJ = Bs[k]
        dFdp += (Bk.T @ f["dsigma_dp"]) * (detJ * wts[k])
    ds_oti = dFdp.T @ w

    def reaction(pvals):
        pr = make_props(contract, pvals); F = np.zeros(ndof)
        for k in range(nip):
            s, _, _ = mat._reg_step(pr, np.zeros(mat.ntens), np.zeros(mat.nstatev),
                                    eps_ip[k], 1.0)
            Bk, detJ = Bs[k]
            F += Bk.T @ s * (detJ * wts[k])
        return F @ w

    rel = []
    for name in names:
        h = 1e-4 * abs(M2_PROPS[name])
        vp = dict(M2_PROPS); vm = dict(M2_PROPS); vp[name] += h; vm[name] -= h
        fd = (reaction(vp) - reaction(vm)) / (2 * h)
        rel.append(abs(ds_oti[names.index(name)] - fd) / (abs(fd) + 1e-30))
    return names, rel


# ===========================================================================
# Test B: full residual solve under force control -- d(tip disp)/dp vs FD
# ===========================================================================
def _newton(mat, props, Xe, etype, fixed_dofs, fext, N=1, dt=1.0):
    pts, wts = k3.gauss_rule(etype); nip = len(wts)
    ndof = 3 * k3.element(etype).nnode
    Bs = [k3.b_matrix(Xe, pts[k], etype) for k in range(nip)]
    free = np.array([d for d in range(ndof) if d not in set(fixed_dofs)], int)
    u = np.zeros(ndof); nt = mat.ntens; ns = mat.nstatev
    sig = np.zeros((nip, nt)); sv = np.zeros((nip, ns)); eps_path = np.zeros((nip, N, nt))
    for inc in range(N):
        up = u.copy(); sp = sig.copy(); svp = sv.copy(); Fext = fext * ((inc + 1) / N)
        for _ in range(40):
            Fint = np.zeros(ndof); K = np.zeros((ndof, ndof))
            for k in range(nip):
                Bk, detJ = Bs[k]
                s, svk, dd = mat._reg_step(props, sp[k], svp[k], Bk @ (u - up), dt)
                sig[k] = s; sv[k] = svk
                Fint += Bk.T @ s * (detJ * wts[k]); K += Bk.T @ dd @ Bk * (detJ * wts[k])
            R = Fint - Fext
            if np.linalg.norm(R[free]) < 1e-9 * (np.linalg.norm(Fext) + 1e-12):
                break
            u[free] += np.linalg.solve(K[np.ix_(free, free)], -R[free])
        for k in range(nip):
            eps_path[k, inc] = Bs[k][0] @ u
    return u, eps_path


def full_solve(mat, contract, etype, pvals, N=1, Fmax=F_TOTAL, distort=0.06,
               load_dir=LOAD_DIR, dt=1.0, hstep=1e-4):
    """Full residual method on a single element under force control:
    assemble K and dR/dp, solve K du/dp = -dR/dp, evaluate d(tip disp)/dp, and
    compare to centered finite differences of the full nonlinear analysis.
    Returns the per-parameter relative error list, or None if the element is a
    reduced-integration single element (rank-deficient / hourglass)."""
    if etype not in FULLY_INTEGRATED:
        return None
    Xe = ref_coords(etype, distort)
    ndof = 3 * k3.element(etype).nnode
    z = Xe[:, 2]; fixed = np.where(z <= z.min() + 1e-6)[0]; loaded = np.where(z >= z.max() - 1e-6)[0]
    fixed_dofs = [3 * n + i for n in fixed for i in range(3)]
    fext = np.zeros(ndof)
    for n in loaded:
        fext[3 * n:3 * n + 3] += (Fmax / len(loaded)) * load_dir
    props = make_props(contract, pvals)
    names = [p["name"] for p in contract["parameters"]]
    resp = lambda u: float(np.mean([u[3 * n:3 * n + 3] @ load_dir for n in loaded]))
    u, eps_path = _newton(mat, props, Xe, etype, fixed_dofs, fext, N, dt)
    pts, wts = k3.gauss_rule(etype); nip = len(wts)
    Bs = [k3.b_matrix(Xe, pts[k], etype) for k in range(nip)]
    free = np.array([d for d in range(ndof) if d not in set(fixed_dofs)], int)
    K = np.zeros((ndof, ndof)); dR = np.zeros((ndof, len(names)))
    for k in range(nip):
        ep = eps_path[k]; dstran = np.vstack([ep[0], np.diff(ep, axis=0)])
        f = replay_fields(mat, props, [list(r) for r in dstran], [dt] * len(dstran))[-1]
        Bk, detJ = Bs[k]
        K += Bk.T @ f["ddsdde"] @ Bk * (detJ * wts[k])
        dR += (Bk.T @ f["dsigma_dp"]) * (detJ * wts[k])
    dudp = np.zeros((ndof, len(names)))
    dudp[free, :] = np.linalg.solve(K[np.ix_(free, free)], -dR[free, :])
    ds_oti = [float(np.mean([dudp[3 * n:3 * n + 3, j] @ load_dir for n in loaded]))
              for j in range(len(names))]
    rel = []
    for j, name in enumerate(names):
        h = hstep * abs(pvals[name])
        vp = dict(pvals); vm = dict(pvals); vp[name] += h; vm[name] -= h
        up, _ = _newton(mat, make_props(contract, vp), Xe, etype, fixed_dofs, fext, N, dt)
        um, _ = _newton(mat, make_props(contract, vm), Xe, etype, fixed_dofs, fext, N, dt)
        fd = (resp(up) - resp(um)) / (2 * h)
        rel.append(abs(ds_oti[j] - fd) / (abs(fd) + 1e-30))
    return rel


# ===========================================================================
def _plot_table(rows, names, fcc=None):
    fig, ax = plt.subplots(figsize=(10.6, 3.5)); ax.axis("off")
    cols = (["Element", "IPs", "DOF", "shape\npatch"]
            + ["assembly\n∂R/∂%s" % n for n in names] + ["full solve\n∂u/∂p"])
    t = ax.table(cellText=[r[1:] for r in rows], colLabels=cols, loc="center", cellLoc="center")
    t.auto_set_font_size(False); t.set_fontsize(10); t.scale(1, 1.6)
    for c in range(len(cols)):
        t[0, c].set_facecolor("#1a2a4a"); t[0, c].set_text_props(color="white", fontweight="bold")
    for ri in range(1, len(rows) + 1):
        t[ri, 0].set_text_props(fontweight="bold")
        for c in range(len(cols)):
            t[ri, c].set_facecolor("#f4f7fb" if ri % 2 else "#e9eef5")
    ax.set_title("Residual method validated across 6 element formulations  (cubic-elastic M2)\n"
                 "assembly: OTI ∂R/∂p vs material FD   ·   full solve: OTI ∂u/∂p vs full-analysis FD",
                 fontsize=12, fontweight="bold", color="#1a2a4a", pad=10)
    if fcc:
        worst = max(v["worst_rel"] for v in fcc.values())
        ets = "/".join(fcc.keys())
        fig.subplots_adjust(bottom=0.16)
        fig.text(0.5, 0.045,
                 "Path-dependent check — FCC crystal plasticity (10 parameters, [100] loading) "
                 "assembled through %s:\nmarch-based DSIGMA_DP vs finite differences, worst error %.1e "
                 "(element-independent, as it must be for a homogeneous state)." % (ets, worst),
                 ha="center", va="bottom", fontsize=9.5, color="#1a5e2a",
                 bbox=dict(boxstyle="round,pad=0.5", fc="#eef7f0", ec="#2e8b57"))
    fig.savefig(os.path.join(FIGDIR, "element_validation_table.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


# FCC crystal-plasticity properties (M6), same as fcc_crystal_results.py
FCC_PROPS = {"g0": 13.0, "h0": 800.0, "q": 1.4, "gd0": 0.001, "m": 0.05,
             "gsat": 55.0, "C11": 168000.0, "C12": 121000.0, "C44": 75000.0, "a": 2.0}
AXIAL = np.array([0.0, 0.0, 1.0])


def fcc_across_elements(N=25, eps_max=0.012, dt=0.04):
    """Path-dependent demonstration: the march-based FCC crystal-plasticity
    DSIGMA_DP (10 parameters, plastic flow) assembled through each NEW element
    type. A homogeneous [100] uniaxial-strain path is prescribed, every IP is
    marched, and the assembled reaction-force sensitivity is compared to centered
    finite differences of the regular UMAT marched over the same path. Inactive
    parameters (C44 has no effect under [100]) use the deck's global-scale
    normalization so their machine-zero sensitivity is not a 0/0 blow-up."""
    md = os.path.join(P1, "oti_provider", "materials", "m6_fcc")
    contract = json.load(open(os.path.join(md, "umat_m6_fcc_oti.json")))
    obj = os.path.join(md, contract["object"]["file"])
    mat = PathMaterial(obj, contract, workdir=os.path.join(FIGDIR, "_m6_link"))
    names = [p["name"] for p in contract["parameters"]]
    props = make_props(contract, FCC_PROPS)
    # approximately uniaxial-STRESS [100] path: axial ramp with lateral contraction
    # (keeps the stress deviatoric / tau-g moderate so the FD reference is well conditioned)
    da = eps_max / N
    de = np.array([-0.35 * da, -0.35 * da, da, 0.0, 0.0, 0.0])
    path = [list(de)] * N; dts = [dt] * N
    print("\n Path-dependent check: FCC plasticity DSIGMA_DP through new elements ([100], %d incr)" % N)
    out = {}
    for et in ["C3D8", "C3D4", "C3D20", "C3D10"]:
        Xe = ref_coords(et, distort=0.0)
        pts, wts = k3.gauss_rule(et); nip = len(wts); ndof = 3 * k3.element(et).nnode
        Bs = [k3.b_matrix(Xe, pts[k], et) for k in range(nip)]
        z = Xe[:, 2]; bottom = np.where(z <= z.min() + 1e-6)[0]
        w = np.zeros(ndof)
        for n in bottom:
            w[3 * n:3 * n + 3] = AXIAL
        # OTI: homogeneous state -> one march, assembled through the element
        dsig = mat.march_fast(props, path, dts)[-1]["dsigma_dp"]  # (6, nparam)
        dFdp = np.zeros((ndof, len(names)))
        for k in range(nip):
            Bk, detJ = Bs[k]
            dFdp += (Bk.T @ dsig) * (detJ * wts[k])
        ds_oti = w @ dFdp

        def reaction(pr):
            s = mat.march_regular(pr, path, dts)[-1]["stress"]
            F = np.zeros(ndof)
            for k in range(nip):
                Bk, detJ = Bs[k]
                F += Bk.T @ s * (detJ * wts[k])
            return F @ w

        ds_fd = np.zeros(len(names))
        for j, nm in enumerate(names):
            h = 1e-4 * abs(FCC_PROPS[nm]) + 1e-12
            vp = dict(FCC_PROPS); vm = dict(FCC_PROPS); vp[nm] += h; vm[nm] -= h
            ds_fd[j] = (reaction(make_props(contract, vp)) - reaction(make_props(contract, vm))) / (2 * h)
        gscale = 1e-3 * np.max(np.abs(ds_fd))                    # deck global-scale normalization
        rel = np.abs(ds_oti - ds_fd) / (np.abs(ds_fd) + gscale)
        out[et] = {"worst_rel": float(np.max(rel)),
                   "rel": {nm: float(r) for nm, r in zip(names, rel)}}
        print("   %-6s IP=%2d  worst ∂R/∂p error over 10 params vs march FD = %.2e"
              % (et, nip, np.max(rel)))
    return out


def main():
    md = os.path.join(P1, "oti_provider", "materials", "m2_cubic")
    contract = json.load(open(os.path.join(md, "umat_m2_cubic_oti.json")))
    obj = os.path.join(md, contract["object"]["file"])
    mat = PathMaterial(obj, contract, workdir=os.path.join(FIGDIR, "_m2_link"))

    print("=" * 78)
    print(" Residual method across element formulations (M2 cubic elastic)")
    print("=" * 78)
    rows = []; jrows = {}; names = None
    for et in ELEMS:
        names, patch_rel = patch_test(mat, contract, et)
        shp = max(k3.verify_shape_functions(et).values())
        solve_rel = full_solve(mat, contract, et, M2_PROPS)
        nip = len(k3.gauss_rule(et)[1]); ndof = 3 * k3.element(et).nnode
        solve_cell = ("%.1e" % max(solve_rel)) if solve_rel is not None else "HG*"
        rows.append([et, et, str(nip), str(ndof), "%.0e" % shp]
                    + ["%.1e" % e for e in patch_rel] + [solve_cell])
        jrows[et] = {"nip": nip, "ndof": ndof, "shape_patch": shp,
                     "assembly_rel": patch_rel, "full_solve_rel": solve_rel}
        print(" %-7s IP=%2d DOF=%3d  patch=%.0e  assembly=%s  solve=%s"
              % (et, nip, ndof, shp, " ".join("%.1e" % e for e in patch_rel), solve_cell))
    worst_a = max(max(jrows[e]["assembly_rel"]) for e in ELEMS)
    solves = [max(jrows[e]["full_solve_rel"]) for e in ELEMS if jrows[e]["full_solve_rel"]]
    print(" worst assembly error = %.2e   worst full-solve error = %.2e" % (worst_a, max(solves)))
    print(" HG* = single reduced-integration element is rank-deficient (hourglass);")
    print("       validated at the assembly level, needs hourglass control / a mesh to solve.")

    fcc = fcc_across_elements()
    _plot_table(rows, names, fcc)
    json.dump({"material": "m2_cubic", "params": names, "elements": jrows,
               "worst_assembly": worst_a, "worst_full_solve": max(solves),
               "fcc_plasticity": fcc,
               "fcc_worst": max(v["worst_rel"] for v in fcc.values())},
              open(os.path.join(FIGDIR, "element_validation.json"), "w"), indent=2)
    print(" figure -> element_validation_table.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
