"""The bounded finite-strain C3D8 with Abaqus's mean dilatation, exactly linearised.

``SolidC3D8FiniteStrain`` keeps full integration as its default (every
published number of examples/finite_strain_c3d8 used it, and its Abaqus check
is a homogeneous field where the two coincide). ``options["integration"] =
"mean_dilatation"`` gives what Abaqus C3D8 does on a non-homogeneous field;
its tangent must be the exact derivative of its force (FD over every DOF of a
distorted, rotated, strongly stretched element) and its force the shared
c3d8_nlgeom one.
"""
from __future__ import annotations

import numpy as np
import pytest

from residual_core.formulations import c3d8_kernel as k
from residual_core.formulations.c3d8_nlgeom import C3D8Nlgeom
from residual_core.formulations.solid_c3d8_finite_strain import SolidC3D8FiniteStrain
from residual_core.materials.base import MaterialBinding
from residual_core.materials.neo_hookean import CompressibleNeoHookean

pytestmark = pytest.mark.unit

RNG = np.random.default_rng(42)
X = k.unit_cube_Xe() + RNG.normal(0.0, 0.07, (8, 3))
ROT = np.array([[np.cos(0.7), -np.sin(0.7), 0.0], [np.sin(0.7), np.cos(0.7), 0.0], [0.0, 0.0, 1.0]])
U = ((X @ (np.array([[1.3, 0.2, -0.1], [0.05, 0.9, 0.12], [0.0, -0.08, 1.1]]).T)) @ ROT.T
     + RNG.normal(0.0, 0.04, (8, 3)) - X).ravel()
BINDING = MaterialBinding(CompressibleNeoHookean(), [2.3, 4.1])


def _eval(u, integration):
    return SolidC3D8FiniteStrain().eval_element(1, "C3D8", X, u, {}, None, BINDING, (0, 1), 1,
                                                None, {"integration": integration})


@pytest.mark.parametrize("integration", ["full", "mean_dilatation"])
def test_the_tangent_is_the_derivative_of_the_force(integration):
    _, K, _, _ = _eval(U, integration)
    errors = []
    for h in (1e-5, 1e-6, 1e-7):
        fd = np.zeros((24, 24))
        for j in range(24):
            e = np.zeros(24)
            e[j] = h
            fd[:, j] = (_eval(U + e, integration)[0] - _eval(U - e, integration)[0]) / (2 * h)
        errors.append(np.max(np.abs(K - fd)) / np.max(np.abs(fd)))
    assert min(errors) < 1e-8, errors


def test_mean_dilatation_force_is_the_shared_abaqus_contract_force():
    r, _, _, info = _eval(U, "mean_dilatation")
    el = C3D8Nlgeom(X[None], "mean_dilatation")
    Fbar = el.Fbar(U.reshape(1, 8, 3))[0]
    sigma = np.array([BINDING.material.evaluate({"F0": np.eye(3), "F1": F, "element": 1, "ip": q + 1},
                                                np.zeros(0), BINDING, (0, 1), 1, None, {})[0]
                      for q, F in enumerate(Fbar)])
    assert np.allclose(r, el.force(U.reshape(1, 8, 3), sigma[None])[0], rtol=1e-13, atol=1e-13)
    assert info["integration"] == "mean_dilatation"


def test_full_integration_remains_the_default_and_differs_here():
    r_default = SolidC3D8FiniteStrain().eval_element(1, "C3D8", X, U, {}, None, BINDING, (0, 1), 1,
                                                     None, {})[0]
    assert np.array_equal(r_default, _eval(U, "full")[0])
    assert np.max(np.abs(r_default - _eval(U, "mean_dilatation")[0])) > 1e-4 * np.max(np.abs(r_default))


# ---------------------------------------------------------------- R1 (B3 review)
# The mean-dilatation path once cast the stress to float: with Dual1 constants
# it returned a real residual (the derivative silently dropped), with OTILib
# ones it raised. Quantity: the element residual r(u; mu, lambda) and tangent
# K, differentiated w.r.t. mu at fixed real u, u_prev and geometry; local (one
# evaluation, no history). Reference: central FD of the real path, 3 steps.

def _fd_mu(h, integration, which):
    plus = MaterialBinding(CompressibleNeoHookean(), [2.3 + h, 4.1])
    minus = MaterialBinding(CompressibleNeoHookean(), [2.3 - h, 4.1])
    run = lambda b: SolidC3D8FiniteStrain().eval_element(
        1, "C3D8", X, U, {}, None, b, (0, 1), 1, None, {"integration": integration})[which]
    return (run(plus) - run(minus)) / (2 * h)


def _fd_check(derivative, integration, which):
    errors = [np.max(np.abs(derivative - _fd_mu(h, integration, which))) for h in (1e-3, 1e-4, 1e-5)]
    scale = np.max(np.abs(_fd_mu(1e-4, integration, which)))
    assert scale > 1e-3                      # a non-trivial derivative
    assert min(errors) < 1e-8 * max(scale, 1.0), errors
    return errors


@pytest.mark.parametrize("integration", ["full", "mean_dilatation"])
def test_a_dual1_parameter_derivative_survives_the_mean_dilatation_path(integration):
    from residual_core.algebra.dual1 import Dual1
    binding = MaterialBinding(CompressibleNeoHookean(), [Dual1(2.3, 1.0), 4.1])
    r, K, _, _ = SolidC3D8FiniteStrain().eval_element(
        1, "C3D8", X, U, {}, None, binding, (0, 1), 1, None, {"integration": integration})
    assert all(isinstance(v, Dual1) for v in np.ravel(r))
    real = _eval(U, integration)
    assert np.allclose([v.real for v in r], real[0], rtol=1e-13, atol=1e-13)
    assert np.allclose([[v.real for v in row] for row in K], real[1], rtol=1e-12, atol=1e-12)
    _fd_check(np.array([v.imag for v in r]), integration, 0)
    _fd_check(np.array([[v.imag for v in row] for row in K]), integration, 1)


@pytest.mark.parametrize("integration", ["full", "mean_dilatation"])
def test_an_oti_parameter_derivative_survives_the_mean_dilatation_path(integration):
    from residual_core.algebra import otilib_adapter
    if not otilib_adapter.otilib_available():
        pytest.skip("OTILib not installed")
    ctx = otilib_adapter.OtiContext(1, 1)
    binding = MaterialBinding(CompressibleNeoHookean(), [ctx.seed(2.3, 1), 4.1])
    r = SolidC3D8FiniteStrain().eval_element(
        1, "C3D8", X, U, {}, None, binding, (0, 1), 1, None,
        {"integration": integration, "compute_tangent": False})[0]
    exponents = ctx.order_directions(1)[0]["exponents"]
    d = np.array([ctx.coeff(v, exponents) for v in r])
    _fd_check(d, integration, 0)
