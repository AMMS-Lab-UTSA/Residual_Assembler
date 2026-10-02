"""The NLGEOM C3D8 operators and their exact linearisation, without a UMAT.

Every derivative the history engine seeds into the material, and the
fixed-stress derivative of the residual, is compared with a centred difference
of the same function (relative step 1e-6; agreement ~1e-10 is truncation +
round-off, a wrong term shows at O(1)). Plus the structural properties the
NLGEOM contract rests on: a pure rotation increment has DSTRAN = 0 and
DROT = that rotation exactly (Hughes-Winget is the Cayley map), DROT is proper
orthogonal, the full-integration force is RA's existing finite-strain force,
the centroid force is the B-bar force of residual_core.replay.kinematics, and
the force of a uniform stress is self-equilibrated.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from residual_core.formulations import c3d8_kernel as k  # noqa: E402
from residual_core.formulations import c3d8_nlgeom as nl  # noqa: E402

pytestmark = pytest.mark.unit

RNG = np.random.default_rng(20261001)
XE = np.array([k.unit_cube_Xe() + RNG.normal(0, 0.06, (8, 3)) for _ in range(3)])
U0 = RNG.normal(0, 0.05, (3, 8, 3))
U1 = U0 + RNG.normal(0, 0.04, (3, 8, 3))
SIG = RNG.normal(0, 1.0, (3, 8, 6))
H = 1e-6


def _rotation(axis, angle):
    from residual_core.corpus.mesh import rotation_matrix
    return rotation_matrix(axis, angle)


def _fd(fun, x, dx):
    return (fun(x + H * dx) - fun(x - H * dx)) / (2 * H)


def _rel(a, b):
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))


@pytest.mark.parametrize("integration", nl.INTEGRATIONS)
def test_the_increment_derivatives_are_exact(integration):
    el = nl.C3D8Nlgeom(XE, integration)
    d0 = RNG.normal(0, 1, (3, 1, 8, 3))
    d1 = RNG.normal(0, 1, (3, 1, 8, 3))
    inc = el.increment(U0, U1)
    d = el.increment_direction(inc, d0, d1)
    plus = el.increment(U0 + H * d0[:, 0], U1 + H * d1[:, 0])
    minus = el.increment(U0 - H * d0[:, 0], U1 - H * d1[:, 0])
    for mine, key in (("dstran", "dstran"), ("ddrot", "drot"), ("dFbar0", "Fbar0"),
                      ("dFbar1", "Fbar1")):
        fd = (plus[key] - minus[key]) / (2 * H)
        assert _rel(d[mine][:, :, 0], fd) < 1e-7, (integration, mine)


@pytest.mark.parametrize("integration", nl.INTEGRATIONS)
def test_the_fixed_stress_derivative_of_the_force_is_exact(integration):
    el = nl.C3D8Nlgeom(XE, integration)
    d = RNG.normal(0, 1, (3, 2, 8, 3))
    analytic = el.dforce_fixed_stress(U1, SIG, d)
    for j in range(2):
        fd = _fd(lambda U: el.force(U, SIG), U1, d[:, j])
        assert _rel(analytic[:, :, j], fd) < 1e-7


@pytest.mark.parametrize("integration", nl.INTEGRATIONS)
def test_a_pure_rotation_increment_strains_nothing_and_rotates_exactly(integration):
    el = nl.C3D8Nlgeom(XE, integration)
    Q = _rotation((1.0, -2.0, 0.5), 0.7)
    x0 = XE + U0
    U_rot = x0 @ Q.T - XE
    inc = el.increment(U0, U_rot)
    assert np.max(np.abs(inc["dstran"])) < 1e-14
    assert np.max(np.abs(inc["drot"] - Q)) < 1e-14
    assert np.max(np.abs(inc["Fbar1"] - Q @ inc["Fbar0"])) < 1e-14


def test_drot_is_proper_orthogonal():
    inc = nl.C3D8Nlgeom(XE).increment(U0, U1)
    R = inc["drot"]
    assert np.max(np.abs(R @ np.swapaxes(R, -1, -2) - np.eye(3))) < 1e-14
    assert np.max(np.abs(np.linalg.det(R) - 1.0)) < 1e-14


def test_full_integration_is_ras_finite_strain_force_and_its_fixed_stress_jacobian():
    el = nl.C3D8Nlgeom(XE, "full")
    ref = np.array([k.element_internal_force_finite_strain(XE[e], U1[e], SIG[e]) for e in range(3)])
    assert _rel(el.force(U1, SIG), ref) < 1e-14
    K = el.dforce_fixed_stress(U1, SIG, el.unit_directions(3))[0]
    assert _rel(K, k.force_tangent_fixed_sigma(XE[0], U1[0], SIG[0])) < 1e-13


@pytest.mark.parametrize("volume,integration", [("mean", "mean_dilatation"), ("centroid", "centroid")])
def test_the_replay_operators_are_this_element(volume, integration):
    from residual_core.replay import kinematics as RK
    el = nl.C3D8Nlgeom(XE, integration)
    B, w = RK.spatial_operators(XE, U1, integration=integration)
    assert _rel(el.force(U1, SIG), np.einsum("eqai,eqa,eq->ei", B, SIG, w)) < 1e-14
    assert _rel(el.Fbar(U1), RK.deformation_gradients(XE, U1, fbar=True, volume=volume)) < 1e-14


def test_the_replay_default_is_mean_dilatation():
    from residual_core.replay import kinematics as RK
    el = nl.C3D8Nlgeom(XE, "mean_dilatation")
    assert _rel(RK.deformation_gradients(XE, U1), el.Fbar(U1)) < 1e-15
    assert _rel(RK.internal_force(XE, U1, SIG), el.force(U1, SIG)) < 1e-15


@pytest.mark.parametrize("integration", nl.INTEGRATIONS)
def test_a_uniform_stress_gives_a_self_equilibrated_element_force(integration):
    el = nl.C3D8Nlgeom(XE, integration)
    f = el.force(U1, np.repeat(SIG[:, :1], 8, axis=1)).reshape(3, 8, 3)
    assert np.max(np.abs(f.sum(axis=1))) < 1e-13 * np.max(np.abs(f))


@pytest.mark.parametrize("integration", nl.INTEGRATIONS)
def test_bbar_matrix_and_force_agree(integration):
    el = nl.C3D8Nlgeom(XE, integration)
    B, w, _ = el.bbar(U1)
    assert _rel(np.einsum("eqsi,eqs,eq->ei", B, SIG, w), el.force(U1, SIG)) < 1e-14
