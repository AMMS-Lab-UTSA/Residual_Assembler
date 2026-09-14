"""dR/dp and du/dp across a mesh, with boundary conditions and a real solve.

The element-level parameter check next door establishes that integrating a
stress sensitivity gives the sensitivity of the integral. It stops there: one
element, no boundary conditions, no scatter, no solve. Everything a parameter
sensitivity is actually used for lives past that point --

    K du/dp = -dR/dp     on the FREE degrees of freedom only,

where the partition decides which rows are solved and which carry a reaction,
the scatter decides where each element's contribution lands, and the sign
decides whether the answer points the right way. A wrong partition produces a
du/dp that is plausible, smooth, and wrong; a sign error produces one that is
exactly wrong and looks like a modelling disagreement.

So this solves the model. Equilibrium is found by Newton at ``p``, the
sensitivity is assembled and solved, and the answer is compared against
re-solving the whole model at ``p + h`` and ``p - h`` and differencing the
converged displacements -- over a plateau of step sizes, because one step
cannot tell a truncation error from a cancellation one.

The ``dsigma/dp`` on the analytic side is DERIVED here from the radial-return
equations, not taken from the same replay that produces the stress, so the two
sides of the comparison are not one implementation checked against itself. The
material is the published SoftwareX J2 reference; what is on trial is this
repository's assembly, partition and solve.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

pytest.importorskip("umat_oti")

from umat_oti.validation.j2_reference import (J2Parameters,  # noqa: E402
                                              J2State, elastic_stiffness,
                                              integrate_increment)

from residual_core.core.finite_difference import sweep_steps  # noqa: E402
from residual_core.formulations.c3d8_kernel import (  # noqa: E402
    ABAQUS_C3D8_GAUSS, b_matrix_reference)
from residual_core.formulations.c3d8_sensitivity import (  # noqa: E402
    apply_dirichlet, assemble_dR_dp, element_dR_dp, solve_du_dp)

#: Relative step sizes on the parameter. Five decades.
STEPS = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6)

#: Measured on this model: the plateau sits near 1e-9, so this fails on a
#: wrong assembly rather than on rounding.
AGREES = 1e-6


# --------------------------------------------------------------------------- #
# the material point: J2 from a frozen incoming state, and its p-derivative
# --------------------------------------------------------------------------- #
def _deviator(sigma):
    pressure = (sigma[0] + sigma[1] + sigma[2]) / 3.0
    return np.array([sigma[0] - pressure, sigma[1] - pressure,
                     sigma[2] - pressure, sigma[3], sigma[4], sigma[5]])


def _mises(dev):
    return float(np.sqrt(1.5 * np.sum(dev[:3] ** 2) + 3.0 * np.sum(dev[3:] ** 2)))


def update(params: J2Parameters, dstrain):
    """One increment from the undeformed state: stress, tangent, whether it
    yielded. The whole strain is the increment, so a displacement field maps
    straight onto a material point without a history to march."""
    result = integrate_increment(params, J2State(), np.asarray(dstrain, float))
    return (np.asarray(result.stress, dtype=float),
            np.asarray(result.ddsdde, dtype=float), result.yielded,
            float(result.statev[0]))


def dsigma_dp_closed_form(params: J2Parameters, dstrain, name: str):
    """``d sigma / d p`` for radial return with linear isotropic hardening.

    From the update, differentiated by hand. With ``q_n = 0`` (one increment
    from rest) and ``dgamma = (q_trial - SIGY0) / (3 mu + H)``::

        d sigma / d SIGY0 = dev_trial * 3 mu / (q_trial (3 mu + H))
        d sigma / d H     = dev_trial * 3 mu dgamma / (q_trial (3 mu + H))

    and in the elastic branch both are exactly zero, because the yield surface
    the parameters describe is not reached.
    """
    C = np.asarray(elastic_stiffness(params), dtype=float)
    trial = C @ np.asarray(dstrain, dtype=float)
    dev_trial = _deviator(trial)
    q_trial = _mises(dev_trial)
    mu = params.E / (2.0 * (1.0 + params.nu))
    if q_trial - params.SIGY0 <= 0.0 or q_trial <= 0.0:
        return np.zeros(6)
    dgamma = (q_trial - params.SIGY0) / (3.0 * mu + params.H)
    common = 3.0 * mu / (q_trial * (3.0 * mu + params.H))
    if name == "SIGY0":
        return dev_trial * common
    if name == "H":
        return dev_trial * common * dgamma
    raise ValueError(f"no closed form here for {name!r}")


# --------------------------------------------------------------------------- #
# a two-element model, assembled by hand so the test owns its own mesh
# --------------------------------------------------------------------------- #
def two_element_mesh():
    coords, index = [], {}
    nid = 1
    for k in (0.0, 1.0):
        for j in (0.0, 1.0):
            for i, x in enumerate((0.0, 1.0, 2.2)):
                coords.append((x, j, k))
                index[(i, int(j), int(k))] = nid
                nid += 1
    connectivity = []
    for e in range(2):
        connectivity.append((e + 1, [
            index[(e, 0, 0)], index[(e + 1, 0, 0)], index[(e + 1, 1, 0)],
            index[(e, 1, 0)], index[(e, 0, 1)], index[(e + 1, 0, 1)],
            index[(e + 1, 1, 1)], index[(e, 1, 1)]]))
    return list(range(1, len(coords) + 1)), np.array(coords, float), connectivity


NODES, COORDS, CONNECTIVITY = two_element_mesh()
NDOF = 3 * len(NODES)
ROW = {nid: 3 * (nid - 1) for nid in NODES}


def element_rows(conn):
    return np.array([ROW[nid] + c for nid in conn for c in range(3)], dtype=int)


def assemble(params: J2Parameters, U: np.ndarray, tangent: bool = True):
    """Global residual (and stiffness) of the two-element model at ``U``."""
    R = np.zeros(NDOF)
    K = np.zeros((NDOF, NDOF)) if tangent else None
    yielded_points = 0
    for eid, conn in CONNECTIVITY:
        Xe = COORDS[[nid - 1 for nid in conn]]
        rows = element_rows(conn)
        ue = U[rows]
        for point, weight in zip(ABAQUS_C3D8_GAUSS.points,
                                 ABAQUS_C3D8_GAUSS.weights):
            B, detJ = b_matrix_reference(Xe, point)
            stress, D, yielded, _q = update(params, B @ ue)
            yielded_points += int(yielded)
            R[rows] += (B.T @ stress) * detJ * weight
            if tangent:
                K[np.ix_(rows, rows)] += (B.T @ D @ B) * detJ * weight
    return R, K, yielded_points


def stress_sensitivity(params: J2Parameters, U: np.ndarray, name: str):
    """``dsigma/dp`` at every integration point of every element."""
    out = {}
    for eid, conn in CONNECTIVITY:
        Xe = COORDS[[nid - 1 for nid in conn]]
        ue = U[element_rows(conn)]
        columns = []
        for point in ABAQUS_C3D8_GAUSS.points:
            B, _detJ = b_matrix_reference(Xe, point)
            columns.append(dsigma_dp_closed_form(params, B @ ue, name))
        out[eid] = np.array(columns)
    return out


def held_and_pulled(stretch: float):
    """x = 0 held in every direction; x = 2.2 pulled along x."""
    prescribed = {}
    for nid in NODES:
        x = COORDS[nid - 1][0]
        if x == 0.0:
            for c in range(3):
                prescribed[ROW[nid] + c] = 0.0
        elif abs(x - 2.2) < 1e-12:
            prescribed[ROW[nid] + 0] = stretch
    return prescribed


def solve_equilibrium(params: J2Parameters, prescribed: dict,
                      tolerance: float = 1e-12):
    """Newton on the free rows, with the prescribed rows held exactly."""
    U = np.zeros(NDOF)
    for dof, value in prescribed.items():
        U[dof] = value
    free = np.array([d for d in range(NDOF) if d not in prescribed], dtype=int)
    for _ in range(40):
        R, K, yielded = assemble(params, U)
        scale = max(float(np.max(np.abs(R))), 1.0)
        if float(np.max(np.abs(R[free]))) / scale < tolerance:
            return U, free, yielded
        U[free] += np.linalg.solve(K[np.ix_(free, free)], -R[free])
    raise AssertionError("the two-element model did not converge")


# --------------------------------------------------------------------------- #
# the checks
# --------------------------------------------------------------------------- #
@pytest.mark.integration
@pytest.mark.parametrize("name", ["SIGY0", "H"])
def test_the_global_parameter_derivative_is_the_derivative_of_the_residual(name):
    """dR/dp assembled and scattered, against a difference of R itself.

    At a FIXED displacement, so what is differenced is the explicit
    dependence -- the one the sensitivity solve needs on its right-hand side.
    """
    params = J2Parameters()
    stretch = 0.01
    U, free, yielded = solve_equilibrium(params, held_and_pulled(stretch))
    assert yielded > 0, (
        "every integration point stayed elastic, so dsigma/dp is identically "
        "zero and this test would pass with the assembly deleted")

    dR_dp, _index = assemble_dR_dp(NODES, COORDS, CONNECTIVITY,
                                   stress_sensitivity(params, U, name))
    assert float(np.max(np.abs(dR_dp))) > 0.0

    reference = float(dict(zip(("E", "NU", "SIGY0", "H"),
                               params.as_tuple()))[name])

    def difference(h):
        plus, _K, _y = assemble(params.with_replaced(name, reference + h), U,
                                tangent=False)
        minus, _K, _y = assemble(params.with_replaced(name, reference - h), U,
                                 tangent=False)
        return (plus - minus) / (2.0 * h)

    sweep = sweep_steps(dR_dp, difference, steps=STEPS, tolerance=AGREES,
                        scale=abs(reference))
    assert not sweep.flat, f"dR/d{name}: {sweep.verdict()}"
    assert sweep.converged, (
        f"dR/d{name} assembled across two elements is not the derivative of "
        f"the assembled residual: {sweep.verdict()}\n"
        f"worst components (global row): {sweep.worst_components()}")


@pytest.mark.integration
@pytest.mark.parametrize("name", ["SIGY0", "H"])
def test_du_dp_from_the_implicit_function_theorem_matches_re_solving(name):
    """The whole sensitivity path, against the only thing that can check it.

    ``K du/dp = -dR/dp`` on the free rows is a claim about a model, and the
    model can be re-solved. Re-solving shares nothing with the sensitivity
    path except the residual itself: not the partition, not the sign, not the
    linear solve.
    """
    params = J2Parameters()
    prescribed = held_and_pulled(0.01)
    U, free, yielded = solve_equilibrium(params, prescribed)
    assert yielded > 0

    _R, K, _y = assemble(params, U)
    dR_dp, _index = assemble_dR_dp(NODES, COORDS, CONNECTIVITY,
                                   stress_sensitivity(params, U, name))
    result = solve_du_dp(K, dR_dp, prescribed.keys())
    assert result.residual_norm < 1e-6 * max(
        float(np.linalg.norm(dR_dp)), 1.0), (
        "the linear solve did not converge; a near-singular reduced stiffness "
        "returns a plausible answer and this is what says it did not")
    assert np.all(result.du_dp[np.array(sorted(prescribed))] == 0.0), (
        "a prescribed displacement does not move when a material parameter "
        "does, so its sensitivity is exactly zero, not nearly")

    reference = float(dict(zip(("E", "NU", "SIGY0", "H"),
                               params.as_tuple()))[name])

    def difference(h):
        plus, _f, _y = solve_equilibrium(
            params.with_replaced(name, reference + h), prescribed)
        minus, _f, _y = solve_equilibrium(
            params.with_replaced(name, reference - h), prescribed)
        return (plus - minus) / (2.0 * h)

    sweep = sweep_steps(result.du_dp, difference, steps=STEPS,
                        tolerance=AGREES, scale=abs(reference))
    assert not sweep.flat, f"du/d{name}: {sweep.verdict()}"
    assert sweep.converged, (
        f"du/d{name} from the sensitivity solve does not match re-solving the "
        f"model: {sweep.verdict()}\n"
        f"worst components (global row): {sweep.worst_components()}")


@pytest.mark.integration
def test_the_sign_of_the_sensitivity_solve_is_the_one_that_is_documented():
    """A sign error passes every symmetry check and is exactly wrong.

    Raising the yield stress of a displacement-driven block cannot increase
    the displacement anywhere -- the displacement is prescribed -- but it
    raises the reaction. The sensitivity of the free rows is what changes, and
    its direction is checkable against a re-solve rather than against taste.
    """
    params = J2Parameters()
    prescribed = held_and_pulled(0.01)
    U, _free, _y = solve_equilibrium(params, prescribed)
    _R, K, _y = assemble(params, U)
    dR_dp, _index = assemble_dR_dp(NODES, COORDS, CONNECTIVITY,
                                   stress_sensitivity(params, U, "SIGY0"))
    du_dp = solve_du_dp(K, dR_dp, prescribed.keys()).du_dp

    step = 1.0
    plus, _f, _y = solve_equilibrium(params.with_replaced("SIGY0", 251.0),
                                     prescribed)
    minus, _f, _y = solve_equilibrium(params.with_replaced("SIGY0", 249.0),
                                      prescribed)
    differenced = (plus - minus) / (2.0 * step)
    moving = np.abs(differenced) > 1e-3 * float(np.max(np.abs(differenced)))
    assert moving.any(), "nothing moved, so there is no sign to check"
    assert np.all(np.sign(du_dp[moving]) == np.sign(differenced[moving])), (
        "the assembled sensitivity points the opposite way from re-solving: "
        "solve_du_dp applies the minus in K du/dp = -dR/dp, and something "
        "here applies it twice or not at all")


@pytest.mark.integration
def test_a_constrained_dof_outside_the_system_is_refused():
    """The partition is checked rather than trusted: an index that is not a
    degree of freedom would otherwise silently eliminate the wrong row."""
    K = np.eye(6)
    with pytest.raises(ValueError) as raised:
        apply_dirichlet(K, np.zeros(6), [0, 9])
    assert "outside the system" in str(raised.value)


@pytest.mark.integration
def test_an_elastic_model_has_no_yield_parameter_sensitivity():
    """The control that says the tests above measure something.

    Below yield the residual does not know SIGY0 exists, so its derivative is
    exactly zero -- not small, zero. A sensitivity path that produced noise
    here would produce noise everywhere and the plateau would hide it.
    """
    params = J2Parameters()
    U, _free, yielded = solve_equilibrium(params, held_and_pulled(1e-5))
    assert yielded == 0, "this stretch was supposed to stay elastic"
    dR_dp, _index = assemble_dR_dp(NODES, COORDS, CONNECTIVITY,
                                   stress_sensitivity(params, U, "SIGY0"))
    assert np.all(dR_dp == 0.0)
