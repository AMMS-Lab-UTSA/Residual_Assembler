"""C3D8 parameter sensitivity: dR/dp, boundary elimination, and du/dp.

The element algebra is already checked by the kernel's own patch tests. What is
checked here is the sensitivity chain built on top of it, against a case whose
answer is known in closed form rather than by finite differences.

For linear isotropic elasticity at fixed Poisson ratio the stiffness is
proportional to Young's modulus, so with a load that does not depend on it

    K(E) u = f,   K = E * Khat   =>   u = (1/E) Khat^-1 f
    du/dE = -u/E   exactly.

That makes every step falsifiable at once: the stress sensitivity fed in, the
B^T integration, the elimination of constrained degrees of freedom, the sign,
and the solve. A mistake in any of them shows up as a mismatch against an exact
result rather than against a tolerance someone chose.
"""

from __future__ import annotations

import numpy as np
import pytest

from residual_core.formulations.c3d8_kernel import (
    ABAQUS_C3D8_GAUSS, b_matrix_reference, element_tangent, isotropic_D,
    unit_cube_Xe,
)
from residual_core.formulations.c3d8_sensitivity import (
    apply_dirichlet, assemble_dR_dp, element_dR_dp, solve_du_dp,
)

E_MODULUS = 200000.0
POISSON = 0.3


def _one_element_system(E=E_MODULUS, nu=POISSON):
    """A single C3D8 unit cube, encastre on z=0, pulled by nodal forces on z=1."""
    Xe = unit_cube_Xe()
    D = isotropic_D(E, nu)
    K = element_tangent(Xe, np.zeros((8, 3)), [D] * 8, mode="small")

    bottom = [a for a in range(8) if Xe[a, 2] < 0.5]
    top = [a for a in range(8) if Xe[a, 2] > 0.5]
    constrained = [3 * a + i for a in bottom for i in range(3)]

    f = np.zeros(24)
    for a in top:
        f[3 * a + 2] = 25.0          # total 100 in +z
    return Xe, D, K, f, constrained, top


def _solve(K, f, constrained):
    Kff, ff, free = apply_dirichlet(K, f, constrained)
    u = np.zeros(24)
    u[free] = np.linalg.solve(Kff, ff)
    return u


def test_element_dR_dp_matches_direct_quadrature():
    """The assembly must be the internal-force integral with dsigma/dp in it."""
    Xe = unit_cube_Xe()
    rng = np.random.default_rng(20260825)
    dsigma = rng.normal(size=(8, 6))
    expected = np.zeros(24)
    for k, weight in enumerate(ABAQUS_C3D8_GAUSS.weights):
        B, detJ = b_matrix_reference(Xe, ABAQUS_C3D8_GAUSS.points[k])
        expected += (B.T @ dsigma[k]) * detJ * weight
    assert np.allclose(element_dR_dp(Xe, dsigma, mode="small"), expected, atol=0, rtol=0)


def test_dR_dp_rejects_the_wrong_shape():
    with pytest.raises(ValueError, match=r"\(8, 6\)"):
        element_dR_dp(unit_cube_Xe(), np.zeros((8, 3)), mode="small")


def test_constrained_dofs_keep_exactly_zero_sensitivity():
    """A prescribed displacement does not move when a material parameter does."""
    Xe, D, K, f, constrained, _top = _one_element_system()
    u = _solve(K, f, constrained)
    sigma_ip = _stress_at_ips(Xe, D, u)
    dR = element_dR_dp(Xe, sigma_ip / E_MODULUS, mode="small")
    result = solve_du_dp(K, dR, constrained)
    assert np.all(result.du_dp[np.asarray(constrained)] == 0.0)


def _stress_at_ips(Xe, D, u):
    sigma = np.zeros((8, 6))
    for k in range(len(ABAQUS_C3D8_GAUSS.weights)):
        B, _detJ = b_matrix_reference(Xe, ABAQUS_C3D8_GAUSS.points[k])
        sigma[k] = D @ (B @ u)
    return sigma


def test_du_dp_reproduces_the_closed_form_for_linear_elasticity():
    """The manufactured case: du/dE = -u/E, exactly.

    sigma = D(E) eps(u) with D proportional to E, so the explicit stress
    sensitivity at fixed displacement is sigma/E, and equilibrium then gives
    du/dE = -u/E with no approximation anywhere.
    """
    Xe, D, K, f, constrained, _top = _one_element_system()
    u = _solve(K, f, constrained)
    assert np.max(np.abs(u)) > 0.0, "the manufactured problem did not deform"

    sigma_ip = _stress_at_ips(Xe, D, u)
    dR_dp = element_dR_dp(Xe, sigma_ip / E_MODULUS, mode="small")
    result = solve_du_dp(K, dR_dp, constrained)

    expected = -u / E_MODULUS
    scale = np.max(np.abs(expected))
    assert scale > 0.0
    assert np.max(np.abs(result.du_dp - expected)) / scale < 1e-10
    assert result.residual_norm < 1e-8 * max(1.0, np.max(np.abs(dR_dp)))


def test_the_explicit_term_equals_the_internal_force_over_E():
    """At equilibrium dR/dE must be f/E, which is what makes the identity work."""
    Xe, D, K, f, constrained, _top = _one_element_system()
    u = _solve(K, f, constrained)
    sigma_ip = _stress_at_ips(Xe, D, u)
    dR_dp = element_dR_dp(Xe, sigma_ip / E_MODULUS, mode="small")
    free = np.setdiff1d(np.arange(24), np.asarray(constrained))
    assert np.allclose(dR_dp[free], f[free] / E_MODULUS, rtol=1e-10, atol=1e-12)


def test_sign_convention_is_documented_and_applied_once():
    """A sign error here is invisible in magnitude and wrong in every use."""
    Xe, D, K, f, constrained, _top = _one_element_system()
    u = _solve(K, f, constrained)
    sigma_ip = _stress_at_ips(Xe, D, u)
    dR_dp = element_dR_dp(Xe, sigma_ip / E_MODULUS, mode="small")
    result = solve_du_dp(K, dR_dp, constrained)
    # Stiffening the material must reduce the displacement it produces.
    top_z = [3 * a + 2 for a in range(8) if Xe[a, 2] > 0.5]
    assert np.all(result.du_dp[top_z] < 0.0)
    assert "sign_convention" in result.as_dict()


def test_global_assembly_matches_the_single_element_case():
    """One element assembled globally must reproduce the element vector."""
    Xe = unit_cube_Xe()
    rng = np.random.default_rng(7)
    dsigma = rng.normal(size=(8, 6))
    node_ids = list(range(1, 9))
    connectivity = [(1, node_ids)]
    assembled, index = assemble_dR_dp(node_ids, Xe, connectivity,
                                      {1: dsigma}, mode="small")
    direct = element_dR_dp(Xe, dsigma, mode="small")
    for local, nid in enumerate(node_ids):
        gi = index[nid]
        assert np.allclose(assembled[3 * gi:3 * gi + 3],
                           direct[3 * local:3 * local + 3])


def test_a_fully_constrained_system_is_reported_not_solved():
    Xe, D, K, f, _constrained, _top = _one_element_system()
    result = solve_du_dp(K, np.ones(24), range(24))
    assert np.all(result.du_dp == 0.0)
    assert "every degree of freedom is constrained" in result.as_dict()["note"]
