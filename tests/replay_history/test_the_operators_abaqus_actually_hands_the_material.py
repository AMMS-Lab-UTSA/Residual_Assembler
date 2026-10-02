"""The finite-strain C3D8 operators, against what Abaqus does.

`batch_operators` builds B once in the REFERENCE configuration, which is
right for a small-strain replay and is why `history_inputs` refuses an
NLGEOM=YES deck outright. These are the operators such a deck needs. They were
first measured on a 160-element Abaqus run (commit 30a5ab4) with the volume
change taken from the element CENTROID -- CENTROID-ERA figures:

    kinematics left at F_ip           stress 26% out, shears right to 0.4%
    kinematics corrected to F_bar     stress 7.3e-05          (centroid F_bar)
    assembly left on plain B          K 28% wrong against a difference of R(u)
    assembly corrected to B_bar       K 3.7e-04               (centroid B_bar)

The second pair is the one worth a test. At the point where K was still 28%
wrong the stress already agreed with the recording to 7e-05 and equilibrium
was already inside Abaqus's own convergence tolerance, so neither of the
checks one would naturally trust said anything.

A later Abaqus probe on a DISTORTED element (B2;
test_mean_dilatation_abaqus_probe.py) showed Abaqus averages the volume change
over the element -- mean dilatation -- and the centroid form is 2.4e-4 off in
DFGRD1 and 2.1e-2 in the reactions there. The two coincide only where J is
affine in the natural coordinates; they are NOT equal in general on regular
bricks either (a random nodal displacement of a 1x0.5x0.5 box separates them at second
order: |Jmean - Jcentroid| 8.6e-5 at amplitude 0.01, 2.5e-3 at 0.05; B3
review). The 7.3e-05 and 3.7e-04 above were therefore measured with the
centroid operators and have NOT been re-measured under mean dilatation: the
cantilever's Abaqus recording and its transformed-UMAT replay are not in this
repository. They stay as the 30a5ab4 measurement, not as a property of the
operators pinned below, which are the mean.
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


def test_the_volume_change_is_the_element_mean():
    """Every point reports the element's mean J, Jbar = sum w0 J / sum w0 =
    v / V -- that is what Abaqus's F_bar means (not the centroid's J)."""
    X = wedge()
    rng = np.random.default_rng(3)
    u = rng.normal(0.0, 0.05, (8, 3))
    bar = kinematics.deformation_gradients(X[None], u[None], fbar=True)[0]
    J, w0 = [], []
    for point, weight in zip(kernel.ABAQUS_C3D8_GAUSS.points, kernel.ABAQUS_C3D8_GAUSS.weights):
        _, detJ0, grad = kernel._b_from_coords(X, point)
        J.append(np.linalg.det(np.eye(3) + u.T @ grad))
        w0.append(detJ0 * weight)
    mean = np.dot(w0, J) / np.sum(w0)
    assert np.allclose([np.linalg.det(F) for F in bar], mean, rtol=1e-11)
    # the 30a5ab4 centroid form is a different number on this element
    _, _, grad_c = kernel._b_from_coords(X, np.zeros(3))
    centroid = np.linalg.det(np.eye(3) + u.T @ grad_c)
    assert abs(centroid - mean) > 1e-6
    bar_c = kinematics.deformation_gradients(X[None], u[None], fbar=True, volume="centroid")[0]
    assert np.allclose([np.linalg.det(F) for F in bar_c], centroid, rtol=1e-11)


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

def test_every_point_shares_the_mean_volumetric_row():
    """The volumetric row is the current-volume average of the points' rows."""
    X = wedge()
    rng = np.random.default_rng(7)
    u = rng.normal(0.0, 0.03, (8, 3))
    B, _ = kinematics.spatial_operators(X[None], u[None])
    rows, vol = [], []
    for point, weight in zip(kernel.ABAQUS_C3D8_GAUSS.points, kernel.ABAQUS_C3D8_GAUSS.weights):
        Bq, detJ = kernel.b_matrix_spatial(X + u, point)
        rows.append(Bq[:3].sum(axis=0))
        vol.append(detJ * weight)
    expected = np.einsum("q,qi->i", vol, rows) / np.sum(vol)
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
    """Mean dilatation at u = 0 IS batch_operators' volume-weighted B-bar:
    the finite-strain operator reduces to the small-strain one exactly."""
    X = wedge()
    zero = np.zeros((1, 8, 3))
    spatial, w_spatial = kinematics.spatial_operators(X[None], zero)
    reference, w_reference = kinematics.batch_operators(X[None])
    assert np.allclose(w_spatial, w_reference, rtol=1e-12)
    assert np.allclose(spatial, reference, atol=1e-12)
    # the centroid form does not reduce to it on a distorted element
    centroid, _ = kinematics.spatial_operators(X[None], zero, integration="centroid")
    assert not np.allclose(centroid[:, :, :3, :], reference[:, :, :3, :], atol=1e-6)


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
