"""dR/dp assembled from DSIGMA_DP, against a difference of R itself.

The bridge carries point-wise parameter sensitivities from the UMAT side and
this repository integrates them:

    dR_e/dp_k = sum_ip  w detJ  B^T (dsigma/dp_k)

The existing end-to-end test checks that sum is arithmetically what it claims
to be. That is worth having and it is not a verification: it compares the
assembly against a hand-written form of the same assembly, so a sign error
shared by both would pass.

This differences the thing itself. The J2 reference is re-run at ``p_k + h``
and ``p_k - h``, the residual is assembled from each of those stress
histories, and the centred difference of the two residuals is compared with
the residual sensitivity assembled from DSIGMA_DP.

What that establishes, exactly: that integrating a sensitivity gives the
sensitivity of the integral -- ``sum B^T dsigma/dp == d/dp sum B^T sigma`` --
through this repository's own quadrature, weights and DOF ordering. It does
NOT establish that ``dsigma/dp`` is right: on this machine the OTI backend is
unavailable, so the pipeline's DSIGMA_DP is itself a difference of the same
reference model, and a wrong dsigma/dp would be wrong on both sides. The UMAT
repository is where dsigma/dp is verified; this is where its assembly is.

Several step sizes, because one cannot tell a truncation error from a
cancellation one -- the same rule the UMAT pipeline holds its own differences
to.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

pytest.importorskip("umat_oti")

from umat_oti.validation.j2_reference import (J2Parameters,  # noqa: E402
                                              build_softwarex_j2_path, run_path)
from umat_oti.validation.parameter_sensitivity import (  # noqa: E402
    ParameterMap, StateMap, compute_j2_parameter_sensitivities)

from residual_core.core.umat_oti_sensitivity import (  # noqa: E402
    integrate_dRe_dp_single)
from residual_core.formulations.c3d8_kernel import (  # noqa: E402
    ABAQUS_C3D8_GAUSS, b_matrix_reference, element_internal_force_small_strain)

UNIT_CUBE = np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
                      [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]])

#: Relative step sizes. Spanning three decades, so a plateau can be told from
#: a coincidence.
STEPS = (1e-3, 1e-4, 1e-5, 1e-6)

#: The parameters the bundled J2 fixture declares, and where each sits in
#: PROPS. This is our own fixture -- ``UMATs/UMATs/generic_ps/j2_props.f`` --
#: so its parameter names and values are published in this project rather than
#: inferred from somebody else's deck.
PARAMETERS = ParameterMap.softwarex_default()
STATES = StateMap.softwarex_default()


def residual_from(stress_history, increment: int) -> np.ndarray:
    """The element internal force at one increment of a stress history."""
    stress = np.tile(np.asarray(stress_history[increment], dtype=float), (8, 1))
    return element_internal_force_small_strain(UNIT_CUBE, stress)


def stresses_at(params: J2Parameters) -> list:
    return [np.asarray(record.stress, dtype=float)
            for record in run_path(params, build_softwarex_j2_path())]


def value_of(params: J2Parameters, name: str) -> float:
    return float(dict(zip(("E", "NU", "SIGY0", "H"), params.as_tuple()))[name.upper()])


def perturbed(params: J2Parameters, name: str, amount: float) -> J2Parameters:
    return params.with_replaced(name, value_of(params, name) + amount)


def assembled_sensitivity(dsigma_dp: np.ndarray) -> np.ndarray:
    """dR_e/dp_k from a point-wise dsigma/dp_k, through the bridge's own sum."""
    B_at_ip, weights, columns = [], [], []
    for point, weight in zip(ABAQUS_C3D8_GAUSS.points, ABAQUS_C3D8_GAUSS.weights):
        B, detJ = b_matrix_reference(UNIT_CUBE, point)
        B_at_ip.append(B)
        weights.append(detJ * weight)
        columns.append(np.asarray(dsigma_dp, dtype=float))
    return integrate_dRe_dp_single(B_at_ip, weights, columns)


@pytest.mark.parametrize("name", ["E", "SIGY0", "H"])
def test_the_assembled_sensitivity_converges_onto_a_difference_of_the_residual(name):
    params = J2Parameters()
    run = compute_j2_parameter_sensitivities(
        params=params, path=build_softwarex_j2_path(),
        parameter_map=PARAMETERS, state_map=STATES)
    index = PARAMETERS.names().index(name)
    # The last increment, which is inside the plastic branch: an elastic one
    # would test a derivative every implementation gets right.
    increment = len(run.increments) - 1
    dsigma_dp = np.asarray(run.increments[increment].dsigma_dp, dtype=float)[:, index]
    assembled = assembled_sensitivity(dsigma_dp)

    reference = value_of(params, name)
    errors = []
    for relative in STEPS:
        step = relative * max(abs(reference), 1.0)
        plus = residual_from(stresses_at(perturbed(params, name, step)), increment)
        minus = residual_from(stresses_at(perturbed(params, name, -step)), increment)
        differenced = (plus - minus) / (2.0 * step)
        scale = max(float(np.max(np.abs(assembled))),
                    float(np.max(np.abs(differenced))), 1e-30)
        errors.append(float(np.max(np.abs(assembled - differenced))) / scale)

    best = min(errors)
    agreeing = [error for error in errors if error < 1e-5]
    assert best < 1e-6, (
        f"dR/d{name} assembled from DSIGMA_DP does not agree with a centred "
        f"difference of the assembled residual: best {best:.3e} over steps "
        f"{STEPS}, errors {['%.2e' % e for e in errors]}")
    assert len(agreeing) >= 2, (
        f"dR/d{name} agreed at only {len(agreeing)} step size(s); one step "
        f"cannot separate truncation error from cancellation. "
        f"errors {['%.2e' % e for e in errors]}")


def test_a_parameter_the_stress_does_not_depend_on_gives_no_residual_change():
    """Poisson's ratio moves the stress, so this is the control that the test
    above is measuring something: a sensitivity that is identically zero must
    assemble to a residual sensitivity that is identically zero."""
    B_at_ip, weights, columns = [], [], []
    for point, weight in zip(ABAQUS_C3D8_GAUSS.points, ABAQUS_C3D8_GAUSS.weights):
        B, detJ = b_matrix_reference(UNIT_CUBE, point)
        B_at_ip.append(B)
        weights.append(detJ * weight)
        columns.append(np.zeros(6))
    assert np.allclose(integrate_dRe_dp_single(B_at_ip, weights, columns), 0.0)


def test_the_state_sensitivity_is_carried_and_is_not_the_stress_one():
    """DSTATEV_DP is a different quantity from DSIGMA_DP and has to survive
    the bridge as one: a run that returned the stress sensitivity twice would
    pass any test that only looked at shapes."""
    run = compute_j2_parameter_sensitivities(
        params=J2Parameters(), path=build_softwarex_j2_path(),
        parameter_map=PARAMETERS, state_map=STATES)
    last = run.increments[-1]
    dstatev = np.asarray(last.dstatev_dp, dtype=float)
    dsigma = np.asarray(last.dsigma_dp, dtype=float)
    assert dstatev.shape == (len(STATES.names()), len(PARAMETERS.names()))
    assert dsigma.shape == (6, len(PARAMETERS.names()))
    # Raising the yield stress cannot increase the accumulated plastic strain.
    assert dstatev[0, PARAMETERS.names().index("SIGY0")] < 0.0
