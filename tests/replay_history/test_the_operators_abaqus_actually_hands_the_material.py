"""The finite-strain C3D8 operators, against what Abaqus does.

`batch_operators` builds B once in the REFERENCE configuration, which is
right for a small-strain replay and is why `history_inputs` refuses an
NLGEOM=YES deck outright. These are the operators such a deck needs, and the
properties pinned here are the ones whose absence was measured on a real
160-element Abaqus run:

    kinematics left at F_ip           stress 26% out, shears right to 0.4%
    kinematics corrected to F_bar     stress 7.3e-05
    assembly left on plain B          K 28% wrong against a difference of R(u)
    assembly corrected to B_bar       K 3.7e-04

The second pair is the one worth a test. At the point where K was still 28%
wrong the stress already agreed with the recording to 7e-05 and equilibrium
was already inside Abaqus's own convergence tolerance, so neither of the
checks one would naturally trust said anything.
"""

import numpy as np
import pytest

from residual_core.formulations import c3d8_kernel as kernel
from residual_core.replay import kinematics


def unit_cube():
    return np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
                     [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]])


def wedge():
    """A trapezoidal element: the centroid and mean-dilatation forms only
    coincide for a parallelepiped, so a distorted element is what separates
    an F_bar that reads the centroid from one that averages."""
    X = unit_cube().copy()
    X[4:, 0] *= 1.6
    X[4:, 2] += 0.35
    return X


def affine(X, gradient):
    return X @ np.asarray(gradient).T


# --------------------------------------------------------------- F_bar

def test_a_uniform_deformation_leaves_the_gradient_alone():
    """F_bar == F wherever J does not vary: nothing to redistribute."""
    X = unit_cube()
    F = np.array([[1.08, 0.04, 0.0], [0.0, 0.95, 0.03], [0.01, 0.0, 1.06]])
    u = affine(X, F - np.eye(3))
    plain = kinematics.deformation_gradients(X[None], u[None], fbar=False)
    bar = kinematics.deformation_gradients(X[None], u[None], fbar=True)
    assert np.allclose(plain, bar, atol=1e-12)
    assert np.allclose(bar[0], F, atol=1e-12)


def test_the_volume_change_comes_from_the_centroid():
    """Every point reports the centroid's J -- that is what F_bar means."""
    X = wedge()
    rng = np.random.default_rng(3)
    u = rng.normal(0.0, 0.05, (8, 3))
    bar = kinematics.deformation_gradients(X[None], u[None], fbar=True)[0]
    _, _, grad_c = kernel._b_from_coords(X, np.zeros(3))
    centroid = np.linalg.det(np.eye(3) + u.T @ grad_c)
    assert np.allclose([np.linalg.det(F) for F in bar], centroid, rtol=1e-11)


def test_the_shape_change_is_untouched():
    """F_bar rescales J and nothing else: the deviatoric part must survive."""
    X = wedge()
    rng = np.random.default_rng(5)
    u = rng.normal(0.0, 0.05, (8, 3))
    plain = kinematics.deformation_gradients(X[None], u[None], fbar=False)[0]
    bar = kinematics.deformation_gradients(X[None], u[None], fbar=True)[0]
    for a, b in zip(plain, bar):
        scale = (np.linalg.det(b) / np.linalg.det(a)) ** (1.0 / 3.0)
        assert np.allclose(b, a * scale, rtol=1e-11)
        # unimodular parts identical
        assert np.allclose(a / np.linalg.det(a) ** (1 / 3),
                           b / np.linalg.det(b) ** (1 / 3), rtol=1e-11)


# --------------------------------------------------------------- B_bar

def test_every_point_shares_the_centroid_volumetric_row():
    X = wedge()
    rng = np.random.default_rng(7)
    u = rng.normal(0.0, 0.03, (8, 3))
    B, _ = kinematics.spatial_operators(X[None], u[None])
    centroid, _ = kernel.b_matrix_spatial(X + u, np.zeros(3))
    expected = centroid[:3].sum(axis=0)
    for q in range(8):
        assert np.allclose(B[0, q, :3, :].sum(axis=0), expected, atol=1e-11)


def test_full_integration_leaves_the_operator_alone():
    X = wedge()
    rng = np.random.default_rng(11)
    u = rng.normal(0.0, 0.03, (8, 3))
    full, w_full = kinematics.spatial_operators(X[None], u[None], integration="full")
    for q, point in enumerate(kernel.ABAQUS_C3D8_GAUSS.points):
        B, detJ = kernel.b_matrix_spatial(X + u, point)
        assert np.allclose(full[0, q], B, atol=1e-12)
        assert np.isclose(w_full[0, q], detJ * kernel.ABAQUS_C3D8_GAUSS.weights[q])


def test_the_shear_rows_are_never_touched():
    """The correction is volumetric. A shear row that moved would be a bug."""
    X = wedge()
    rng = np.random.default_rng(13)
    u = rng.normal(0.0, 0.03, (8, 3))
    full, _ = kinematics.spatial_operators(X[None], u[None], integration="full")
    bar, _ = kinematics.spatial_operators(X[None], u[None])
    assert np.allclose(full[:, :, 3:, :], bar[:, :, 3:, :], atol=1e-14)
    assert not np.allclose(full[:, :, :3, :], bar[:, :, :3, :])


def test_at_no_displacement_it_is_the_reference_operator():
    X = wedge()
    zero = np.zeros((1, 8, 3))
    spatial, w_spatial = kinematics.spatial_operators(X[None], zero)
    reference, w_reference = kinematics.batch_operators(X[None])
    assert np.allclose(w_spatial, w_reference, rtol=1e-12)
    # The two differ only in how the volumetric row is averaged: the centroid
    # here, the volume-weighted mean in batch_operators. They agree on a
    # parallelepiped and differ slightly on this wedge, which is the point.
    assert np.allclose(spatial[:, :, 3:, :], reference[:, :, 3:, :], atol=1e-12)


def test_an_inverted_element_is_refused_rather_than_integrated():
    X = unit_cube()
    u = np.zeros((8, 3))
    u[:, 2] = -2.0 * X[:, 2]                     # turn the element inside out
    with pytest.raises(ValueError, match="inverted|positive Jacobian"):
        kinematics.spatial_operators(X[None], u[None])


# --------------------------------------------------------------- K_geo

def test_the_initial_stress_term_is_symmetric():
    X = wedge()
    rng = np.random.default_rng(17)
    u = rng.normal(0.0, 0.02, (8, 3))
    stress = rng.normal(0.0, 100.0, (1, 8, 6))
    K = kinematics.geometric_stiffness(X[None], u[None], stress)
    assert np.allclose(K[0], K[0].T, atol=1e-9)


def test_no_stress_means_no_initial_stress_term():
    X = wedge()
    u = np.zeros((1, 8, 3))
    K = kinematics.geometric_stiffness(X[None], u, np.zeros((1, 8, 6)))
    assert np.allclose(K, 0.0, atol=1e-14)
