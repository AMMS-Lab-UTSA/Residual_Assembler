#!/usr/bin/env python3
"""Program-2 residual-method sensitivities for the crystal-plasticity flow model
on a 4x4x4 C3D8 mesh under simple shear -- the mesh-independence demonstration.

The residual assembler links the distributed OTI .obj, path-marches the material
at every integration point (carrying DSIGMA_DP/DSTATEV_DP), assembles the global
tangent K and residual sensitivity dR/dp = sum B^T (dsigma/dp) dV, and solves
K (du/dp) = -dR/dp for the free degrees of freedom. Under homogeneous simple
shear the converged field is affine, so du/dp is ~0 and the volume-averaged
sigma_vM sensitivities equal the single-element result: the sensitivities are
mesh independent, exactly as reported in the study.

    UMAT_OTI_ROOT=~/Documents/UMAT_source_transformation \
      python results/cp_mesh_residual.py
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
sys.path.insert(0, RA)
P1 = os.path.expanduser(os.environ.get("UMAT_OTI_ROOT", "~/Documents/UMAT_source_transformation"))
from residual_core.core.model import Element, Model
from residual_core.core.dof_manager import DofManager
from residual_core.formulations import c3d8_kernel as kern
from residual_core.replay.path_material import PathMaterial, mises, dmises_dstress

FIGDIR = os.path.join(RA, "results", "figures"); os.makedirs(FIGDIR, exist_ok=True)
MATDIR = os.path.join(P1, "oti_provider", "materials", "m5_cpflow")
LABELS = [r"$\tau_0$", r"$\Delta G$", r"$q$", r"$p$", r"$\dot{\gamma}_0$", r"$H$"]
COLORS = ["#8e5bd6", "#2a7fd4", "#d64b4b", "#f2c230", "#f08a3c", "#2fa6a0"]


def grid_mesh(n=4, L=1.0):
    """n x n x n C3D8 cube on [0,L]^3. Returns Model + node coord array."""
    npc = n + 1
    nid = lambda i, j, k: 1 + i + npc * (j + npc * k)
    nodes = {}
    for k in range(npc):
        for j in range(npc):
            for i in range(npc):
                nodes[nid(i, j, k)] = (L * i / n, L * j / n, L * k / n)
    elems = {}; eid = 0
    for k in range(n):
        for j in range(n):
            for i in range(n):
                eid += 1
                elems[eid] = [nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                              nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1), nid(i, j+1, k+1)]
    m = Model(nodes={a: tuple(map(float, b)) for a, b in nodes.items()},
              elements={e: Element(e, "C3D8", list(c)) for e, c in elems.items()})
    return m, nodes


def main():
    contract = json.load(open(os.path.join(MATDIR, "umat_m5_cpflow_oti.json")))
    obj = os.path.join(MATDIR, contract["object"]["file"])
    props0 = list(contract["validation"]["props_values"])
    pidx = [p["props_index"] - 1 for p in contract["parameters"]]
    pval = np.array([props0[i] for i in pidx])
    mat = PathMaterial(obj, contract); nprm = mat.nparam
    n = 4; L = 1.0; N = 50; e12_max = 0.05
    m, nodes = grid_mesh(n, L)
    dm = DofManager(m.nodes.keys())
    pts = kern.ABAQUS_C3D8_GAUSS.points; wts = kern.ABAQUS_C3D8_GAUSS.weights

    # Homogeneous simple shear u_x = gamma12 * Y imposed on the WHOLE boundary
    # (prescribed displacement, du/dp = 0); only interior nodes are free. This is
    # the standard controlled-deformation RVE test: the equilibrium field is
    # affine, so the residual solve must return du/dp ~ 0.
    def on_boundary(c):
        return (abs(c[0]) < 1e-9 or abs(c[0] - L) < 1e-9 or abs(c[1]) < 1e-9 or
                abs(c[1] - L) < 1e-9 or abs(c[2]) < 1e-9 or abs(c[2] - L) < 1e-9)
    fixed = set()
    for nid, c in nodes.items():
        if on_boundary(c):
            for d in dm.node_dofs(nid): fixed.add(d)
    free = np.array([d for d in range(dm.ndof) if d not in fixed], int)

    dgam = (2.0 * e12_max) / N
    e12 = np.array([(s + 1) * dgam / 2.0 for s in range(N)])
    dstran = [[0.0, 0.0, 0.0, dgam, 0.0, 0.0] for _ in range(N)]
    dt = [1.0] * N

    # every IP sees the same uniform shear path -> march one representative IP.
    oti = mat.march_oti(props0, dstran, dt)

    # assemble K and dR/dp at each increment; solve K du/dp = -dR/dp (free dofs).
    dudp_norms = []; svm_sens = np.zeros((N, nprm)); svm_curve = np.zeros(N)
    ndof = dm.ndof
    edofs = {el.eid: np.asarray(dm.element_dofs(el.connectivity, ("UX", "UY", "UZ")), int)
             for el in m.elements.values()}
    Bcache = {}
    for el in m.elements.values():
        Xe = m.coords_of(el.connectivity)
        Bcache[el.eid] = [kern.b_matrix_reference(Xe, pts[k]) for k in range(8)]

    for s in range(N):
        D_ip = oti[s]["ddsdde"]; sig = oti[s]["stress"]; dsig_dp = oti[s]["dsigma_dp"]
        K = np.zeros((ndof, ndof)); dR = np.zeros((ndof, nprm))
        for el in m.elements.values():
            ed = edofs[el.eid]; Xe = m.coords_of(el.connectivity)
            Dstack = np.repeat(D_ip[None, :, :], 8, axis=0)
            Ke = kern.element_tangent(Xe=Xe, Ue=None, Dmat_ip=Dstack, mode="small")
            K[np.ix_(ed, ed)] += Ke
            for k in range(8):
                B, detJ = Bcache[el.eid][k]
                dR[ed, :] += (B.T @ dsig_dp) * (detJ * wts[k])
        Kff = K[np.ix_(free, free)]
        dudp = np.zeros((ndof, nprm))
        dudp[free, :] = np.linalg.solve(Kff, -dR[free, :])
        dudp_norms.append(float(np.linalg.norm(dudp)))
        # volume-averaged sigma_vM sensitivity: dsigma/dp per IP = material dsig_dp
        # + D B (du_e/dp). Homogeneous -> du/dp ~ 0, so this equals the material value.
        g = dmises_dstress(sig); vol = 0.0; num = np.zeros(nprm)
        for el in m.elements.values():
            ed = edofs[el.eid]
            for k in range(8):
                B, detJ = Bcache[el.eid][k]
                dsig_ip = dsig_dp + D_ip @ (B @ dudp[ed, :])
                num += g @ dsig_ip * (detJ * wts[k]); vol += detJ * wts[k]
        svm_sens[s] = num / vol
        svm_curve[s] = mises(sig)

    # single-element reference (material point) for the mesh-independence check
    single = np.array([dmises_dstress(o["stress"]) @ o["dsigma_dp"] for o in oti])
    max_mesh_diff = float(np.max(np.abs(svm_sens - single)) / max(np.max(np.abs(single)), 1e-30))

    weighted = np.abs(svm_sens) * np.abs(pval)[None, :]
    wpct = 100.0 * weighted / np.maximum(weighted.sum(1, keepdims=True), 1e-30)
    _plot(e12, svm_curve, wpct, n)

    summary = {"mesh": "%dx%dx%d C3D8" % (n, n, n), "n_elements": len(m.elements),
               "n_dof": dm.ndof, "n_increments": N,
               "max_dudp_norm": float(max(dudp_norms)),
               "mesh_vs_single_element_rel": max_mesh_diff,
               "sigma_vM_final_MPa": float(svm_curve[-1])}
    json.dump(summary, open(os.path.join(FIGDIR, "cp_mesh_summary.json"), "w"), indent=2)
    print("=" * 68)
    print(" CP residual method on %dx%dx%d C3D8 mesh (%d elems, %d dof)"
          % (n, n, n, len(m.elements), dm.ndof))
    print("=" * 68)
    print(" max ||du/dp|| over path          = %.2e  (affine field -> ~0)" % max(dudp_norms))
    print(" mesh vs single-element sens (rel)= %.2e  (mesh-independent)" % max_mesh_diff)
    print(" sigma_vM(E12=%.3f)               = %.1f MPa" % (e12[-1], svm_curve[-1]))
    print(" figure -> %s/weighted_sigvm_sensitivities_4x4x4.png" % FIGDIR)
    return 0 if max_mesh_diff < 1e-6 else 1


def _plot(e12, svm, wpct, n):
    fig, axL = plt.subplots(figsize=(7.2, 4.6)); axR = axL.twinx()
    cum = np.zeros_like(e12)
    for k in range(wpct.shape[1]):
        axR.fill_between(e12, cum, cum + wpct[:, k], color=COLORS[k], alpha=0.85, linewidth=0, label=LABELS[k])
        cum = cum + wpct[:, k]
    axL.plot(e12, svm / 1000.0, color="#7a1f1f", lw=2.4, zorder=5)
    axL.set_xlabel(r"Shear strain $E_{12}$"); axL.set_ylabel(r"Stress $\sigma_{vM}$ (GPa)")
    axR.set_ylabel("Weighted sensitivities (%)"); axR.set_ylim(0, 100); axL.set_xlim(e12[0], e12[-1])
    axL.set_title(r"Weighted $\sigma_{vM}$ Parameter Sensitivities (%d$\times$%d$\times$%d Mesh)" % (n, n, n))
    axR.legend(ncol=3, fontsize=8, loc="lower center", framealpha=0.9)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "weighted_sigvm_sensitivities_4x4x4.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
