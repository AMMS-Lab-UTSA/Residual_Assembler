# -*- coding: utf-8 -*-
"""Residual-method parameter-sensitivity engine on one elastic C3D8.

This is the engine every later material plugs into. Given a converged state it:

  1. assembles the residual's parameter derivative  R^(1) = d r / da
        = integral_Omega  B^T (d sigma / da)  dV
  2. builds the tangent  K = d r / d u
  3. solves the sensitivity system  K U^(1) = -R^(1)   ->   du/da
  4. chain-rules to user outputs:  d g / da = (dg/du) . (du/da) + dg/da

The ONLY material-specific input is  d sigma / da  (the integration-point stress
sensitivity). For isotropic elasticity the stiffness D is linear in E, so
`d sigma / dE = sigma / E` in closed form, and that is what this module uses. For a
real UMAT the same field comes from the OTIS-transformed UMAT (written into SDVs);
see ../../docs/otis_umat_connection.md. The engine code does not change — only the
source of `d sigma / da` changes.

Built on the verified C3D8 kernel (`residual_core/formulations/c3d8_kernel.py`).
Run directly for a self-check against central finite differences of the solve;
`verify_against_abaqus.py` checks the same engine against a live Abaqus solve.
"""
from __future__ import annotations

import numpy as np

from residual_core.formulations import c3d8_kernel as k

G = k.ABAQUS_C3D8_GAUSS
NU = 0.3

# unit-cube C3D8, standard node order (0-based nodes 0..7)
XE = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
               [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)


def dof(node, comp):
    """Global DOF index (node-major): node in 0..7, comp in 0..2."""
    return 3 * node + comp


def _b_all(Xe):
    return [k.b_matrix_reference(Xe, xi)[0] for xi in G.points]   # 8 x (6,24)


def stresses(Xe, D, u):
    """Integration-point stresses (8,6) for displacement field u (24,)."""
    return np.array([D @ (Bk @ u) for Bk in _b_all(Xe)])


def solve(E, fixed, fext):
    """Linear-elastic element solve.

    fixed : {dof_index: prescribed_value}   (symmetry + any prescribed displacement)
    fext  : (24,) external nodal load vector
    returns (u, sigma_ip, K, free_dofs)
    """
    D = k.isotropic_D(E, NU)
    K = k.element_tangent(Xe=XE, Ue=None, Dmat_ip=D, mode="small")
    u = np.zeros(24)
    fix = np.array(sorted(fixed), int)
    for d, v in fixed.items():
        u[d] = v
    free = np.array([d for d in range(24) if d not in fixed], int)
    u[free] = np.linalg.solve(K[np.ix_(free, free)],
                              fext[free] - K[np.ix_(free, fix)] @ u[fix])
    return u, stresses(XE, D, u), K, free


def sensitivity_dE(E, fixed, fext):
    """Residual-method sensitivity of the solution and outputs w.r.t. E.

    Returns a dict with:
      u              converged displacement (24,)
      sig            integration-point stress (8,6)
      dudE           du/dE (24,)               <- from K U^(1) = -R^(1)
      dsig_dE_total  total d sigma / dE per IP (8,6)  <- chain rule
    """
    u, sig, K, free = solve(E, fixed, fext)
    # material-specific ingredient: d sigma / dE at fixed strain (elastic: sigma/E)
    dsig_dE = sig / E
    # R^(1) = d f_int / dE at fixed u = integral B^T (d sigma / dE) dV
    R1 = k.element_internal_force_small_strain(XE, dsig_dE)
    # solve K U^(1) = -R^(1); prescribed-DOF sensitivity is 0 (BC values fixed in E)
    U1 = np.zeros(24)
    U1[free] = np.linalg.solve(K[np.ix_(free, free)], -R1[free])
    # chain rule: total d sigma / dE = d sigma/dE|_u + D B (du/dE)
    D = k.isotropic_D(E, NU)
    B = _b_all(XE)
    dsig_tot = np.array([dsig_dE[i] + D @ (B[i] @ U1) for i in range(8)])
    return dict(u=u, sig=sig, dudE=U1, dsig_dE_total=dsig_tot)


def von_mises(s):
    s11, s22, s33, s12, s13, s23 = s
    return np.sqrt(0.5 * ((s11 - s22) ** 2 + (s22 - s33) ** 2 + (s33 - s11) ** 2)
                   + 3 * (s12 ** 2 + s13 ** 2 + s23 ** 2))


# --------------------------------------------------------------------------- #
# standard boundary sets for the 1/8-symmetry unit cube
# --------------------------------------------------------------------------- #
X0, X1 = [0, 3, 4, 7], [1, 2, 5, 6]
Y0, Z0 = [0, 1, 4, 5], [0, 1, 2, 3]


def symmetry_bcs():
    s = {}
    for n in X0:
        s[dof(n, 0)] = 0.0
    for n in Y0:
        s[dof(n, 1)] = 0.0
    for n in Z0:
        s[dof(n, 2)] = 0.0
    return s


def _fd_dE(E, fixed, fext, h_rel=1e-4):
    h = E * h_rel
    up, sp, _, _ = solve(E + h, fixed, fext)
    um, sm, _, _ = solve(E - h, fixed, fext)
    return (up - um) / (2 * h), (sp - sm) / (2 * h)


def _demo():
    E0 = 210000.0
    # Case A: load control (total Fx = 100 on the x=1 face) -> nonzero du/dE
    fixed = symmetry_bcs()
    fext = np.zeros(24)
    for n in X1:
        fext[dof(n, 0)] = 100.0 / 4
    res = sensitivity_dE(E0, fixed, fext)
    du_fd, _ = _fd_dE(E0, fixed, fext)
    d = dof(1, 0)
    print("CASE A  load control")
    print("  du/dE (loaded node x): engine %+.6e  FD %+.6e  analytic(-u/E) %+.6e"
          % (res["dudE"][d], du_fd[d], -res["u"][d] / E0))

    # Case B: displacement control (u_x=1e-3 on x=1 face) -> nonzero d sigma11/dE
    fixed = symmetry_bcs()
    for n in X1:
        fixed[dof(n, 0)] = 1.0e-3
    res = sensitivity_dE(E0, fixed, np.zeros(24))
    _, ds_fd = _fd_dE(E0, fixed, np.zeros(24))
    print("CASE B  displacement control")
    print("  dSigma11/dE (IP0): engine %+.6e  FD %+.6e  analytic(strain) %+.6e"
          % (res["dsig_dE_total"][0][0], ds_fd[0][0], 1.0e-3))


if __name__ == "__main__":
    _demo()
