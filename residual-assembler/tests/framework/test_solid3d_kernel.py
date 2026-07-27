"""solid3d_kernel: 3D solid element formulations for the residual assembler.

Pins the guarantees the residual method relies on for every supported element
type (C3D8, C3D8R, C3D20, C3D20R, C3D4, C3D10):

  * partition of unity          sum_a N_a = 1
  * linear completeness         sum_a N_a X_a = x(xi)  (reproduces coordinates)
  * constant-gradient           sum_a dN_a/dxi = 0
  * gradient correctness        analytic dN == complex-step dN
  * constant-strain patch test  a linear displacement field on a DISTORTED
                                 element yields exactly the constant strain at
                                 every integration point (the property that makes
                                 the assembled residual sensitivity correct for
                                 an arbitrary element shape)

Also checks a rigid-body translation produces zero strain, and that the tangent
of a linear-elastic material equals a finite difference of the internal force
(K = dFint/du) on a distorted element -- i.e. the assembled B/detJ/weights are
mutually consistent.

Pure numpy; no OTILib, no Abaqus. Run: python tests/framework/test_solid3d_kernel.py
"""
import os
import sys

import numpy as np

RA = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.insert(0, RA)
from residual_core.formulations import solid3d_kernel as k3

TOL = 1e-11


def _fail(msg):
    print("  FAIL:", msg)
    return False


def test_shape_functions():
    ok = True
    for et in k3.SUPPORTED:
        r = k3.verify_shape_functions(et)
        worst = max(r.values())
        tag = "ok" if worst < TOL else "FAIL"
        print("  %-7s pou=%.0e lin=%.0e grad_sum=%.0e grad_num=%.0e patch=%.0e  [%s]"
              % (et, r["pou"], r["lin"], r["grad_sum"], r["grad_num"], r["patch"], tag))
        if worst >= TOL:
            ok = _fail("%s shape-function residual %.2e >= %.0e" % (et, worst, TOL))
    return ok


def test_rigid_body_zero_strain():
    """A pure translation must produce zero strain at every IP (distorted elem)."""
    ok = True
    rng = np.random.RandomState(1)
    for et in k3.SUPPORTED:
        e = k3.element(et)
        nat = {"C3D8": k3._HEX_CORNERS, "C3D8R": k3._HEX_CORNERS, "C3D20": k3.HEX20_NODES,
               "C3D20R": k3.HEX20_NODES, "C3D4": k3.TET4_NODES, "C3D10": k3.TET10_NODES}[et]
        Xe = nat + 0.1 * rng.randn(*nat.shape)
        Ue = np.tile(np.array([0.3, -0.2, 0.7]), (e.nnode, 1)).reshape(-1)   # translation
        worst = 0.0
        for xi in e.points:
            B, _ = k3.b_matrix(Xe, xi, et)
            worst = max(worst, np.max(np.abs(B @ Ue)))
        if worst >= TOL:
            ok = _fail("%s rigid translation strain %.2e" % (et, worst))
    print("  rigid-body translation -> zero strain: %s" % ("ok" if ok else "FAIL"))
    return ok


def _elastic_C(E=210000.0, nu=0.3):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu)); mu = E / (2 * (1 + nu))
    C = np.zeros((6, 6))
    C[:3, :3] = lam
    for i in range(3):
        C[i, i] += 2 * mu
    for i in range(3, 6):
        C[i, i] = mu                    # engineering shear -> mu (not 2 mu)
    return C


def test_tangent_matches_fd():
    """K = dFint/du for a linear-elastic material on a distorted element."""
    ok = True
    C = _elastic_C()
    rng = np.random.RandomState(2)
    for et in k3.SUPPORTED:
        e = k3.element(et); ndof = 3 * e.nnode
        nat = {"C3D8": k3._HEX_CORNERS, "C3D8R": k3._HEX_CORNERS, "C3D20": k3.HEX20_NODES,
               "C3D20R": k3.HEX20_NODES, "C3D4": k3.TET4_NODES, "C3D10": k3.TET10_NODES}[et]
        Xe = nat + 0.05 * rng.randn(*nat.shape)
        pts, wts = k3.gauss_rule(et)

        def fint(u):
            sig = np.array([k3.b_matrix(Xe, pts[k], et)[0] @ u for k in range(len(wts))])
            sig = sig @ C.T
            return k3.element_internal_force(Xe, sig, et)

        u0 = 1e-3 * rng.randn(ndof)
        K = k3.element_tangent(Xe, [C] * len(wts), et)
        h = 1e-6; Kfd = np.zeros((ndof, ndof))
        for j in range(ndof):
            up = u0.copy(); um = u0.copy(); up[j] += h; um[j] -= h
            Kfd[:, j] = (fint(up) - fint(um)) / (2 * h)
        rel = np.max(np.abs(K - Kfd)) / (np.max(np.abs(K)) + 1e-30)
        if rel >= 1e-6:
            ok = _fail("%s tangent vs FD rel %.2e" % (et, rel))
        else:
            print("  %-7s K == dFint/du (rel %.1e)  ok" % (et, rel))
    return ok


def main():
    print("solid3d_kernel element formulations")
    results = [test_shape_functions(), test_rigid_body_zero_strain(), test_tangent_matches_fd()]
    if all(results):
        print("PASS: all element formulations verified")
        return 0
    print("FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
