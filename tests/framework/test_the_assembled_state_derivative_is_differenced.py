"""dR/dq, assembled from a point-wise state sensitivity, against a difference.

The third derivative of the residual, and the one this repository did not
have. ``R(u, p, q)`` depends on the state the material carried INTO the
increment as independently as it depends on the displacement: freeze ``u``,
move the accumulated plastic strain at one integration point, and the
residual moves.

    dR_e/dq_(k,j) = w_k detJ_k B_k^T (d sigma_k / d q_(k,j))

Two things are on trial and they are different questions.

**The integral.** Does integrating a state sensitivity give the sensitivity of
the integral, through this repository's own quadrature, weights, B matrices
and DOF ordering? Answered by differencing the assembled residual with respect
to the state, over a plateau of step sizes.

**The structure.** A parameter is shared by every integration point; a state
is not. ``q_(k,j)`` moves the stress at point ``k`` and nowhere else, so
``dR/dq`` is block-sparse by construction and an assembly that treats the
state as global produces a dense matrix whose diagonal blocks are right and
whose other entries are all wrong -- which a one-column check cannot see.
Answered by asserting the zeros, and by checking that a perturbation at one
element's point moves no row belonging to a node that element does not touch.

The ``d sigma / d q`` differenced against is DERIVED here in closed form from
the radial-return equations rather than taken from the same model that
produces the stress, so the two sides are not one implementation checked
against itself. What this does NOT establish is that a transformed UMAT's own
DSTATEV/DSIGMA sensitivities are right -- that is verified in the UMAT
repository, upstream of this.
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
                                              J2State, build_softwarex_j2_path,
                                              elastic_stiffness,
                                              integrate_increment)

from residual_core.core.finite_difference import sweep_steps  # noqa: E402
from residual_core.core.state_sensitivity import (  # noqa: E402
    StateSensitivityError, assemble_dR_dq, element_dR_dq, state_column)
from residual_core.core.where_it_went_wrong import (  # noqa: E402
    StateSensitivityProbe, diagnose)
from residual_core.formulations.c3d8_kernel import (  # noqa: E402
    ABAQUS_C3D8_GAUSS, b_matrix_reference, element_internal_force_small_strain)
from verified_fixtures import UNIT_CUBE, hex_fixtures  # noqa: E402

#: Step sizes relative to the state's own size. Six decades, so truncation
#: dominates one end and cancellation the other.
STEPS = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7)

#: Measured: the plateau sits at 1e-9 to 1e-11 for the analytic dsigma/dq
#: below, so this leaves two decades and fails on a wrong assembly.
AGREES = 1e-7


# --------------------------------------------------------------------------- #
# the material point, and its state derivative in closed form
# --------------------------------------------------------------------------- #
def _deviator(sigma):
    pressure = (sigma[0] + sigma[1] + sigma[2]) / 3.0
    return np.array([sigma[0] - pressure, sigma[1] - pressure,
                     sigma[2] - pressure, sigma[3], sigma[4], sigma[5]])


def _mises(dev):
    return float(np.sqrt(1.5 * np.sum(dev[:3] ** 2) + 3.0 * np.sum(dev[3:] ** 2)))


def dsigma_dq_closed_form(params: J2Parameters, state: J2State, dstrain):
    """``d sigma_(n+1) / d q_n`` for radial return with linear hardening.

    Derived from the update rather than differenced from it. With
    ``dgamma = (q_trial - SIGY0 - H q_n) / (3 mu + H)`` and
    ``sigma = dev_trial (1 - 3 mu dgamma / q_trial) + p_trial``, only
    ``dgamma`` depends on ``q_n``::

        d dgamma / d q_n = -H / (3 mu + H)
        d sigma  / d q_n = dev_trial * 3 mu H / (q_trial (3 mu + H))

    and in the elastic branch the state does not enter the stress at all, so
    the derivative is exactly zero -- which is a fact worth testing, not a
    case to skip.
    """
    C = np.asarray(elastic_stiffness(params), dtype=float)
    trial = np.asarray(state.stress, dtype=float) + C @ np.asarray(dstrain,
                                                                   dtype=float)
    dev_trial = _deviator(trial)
    q_trial = _mises(dev_trial)
    mu = params.E / (2.0 * (1.0 + params.nu))
    yield_stress = params.SIGY0 + params.H * state.statev[0]
    if q_trial - yield_stress <= 0.0 or q_trial <= 0.0:
        return np.zeros(6), False
    factor = 3.0 * mu * params.H / (q_trial * (3.0 * mu + params.H))
    return dev_trial * factor, True


def a_yielded_state(params: J2Parameters = None):
    """Walk the published J2 path until the point is plastic, and stop there.

    A state derivative taken in the elastic branch is zero, and a test that
    only ever looked at zero would pass with the assembly deleted.
    """
    params = params or J2Parameters()
    path = build_softwarex_j2_path()
    state = J2State()
    increments = list(getattr(path, "increments", path))
    for index, dstrain in enumerate(increments):
        result = integrate_increment(params, state, dstrain)
        state = J2State(stress=result.stress, statev=result.statev,
                        stran=result.stran)
        if result.yielded and index + 1 < len(increments):
            return params, state, np.asarray(increments[index + 1], dtype=float)
    pytest.skip("the published J2 path never yields, so there is no plastic "
                "state to take a state derivative at")


def stress_at(params, state, dstrain, eqplas):
    """The stress this increment produces from a given incoming state."""
    moved = J2State(stress=state.stress, statev=(float(eqplas),),
                    stran=state.stran)
    return np.asarray(integrate_increment(params, moved, dstrain).stress,
                      dtype=float)


# --------------------------------------------------------------------------- #
# one element
# --------------------------------------------------------------------------- #
def test_the_closed_form_state_derivative_is_the_derivative_of_the_stress():
    """Before assembling it, check the thing being assembled.

    A wrong ``dsigma/dq`` integrated perfectly would produce a ``dR/dq`` that
    differences perfectly against a residual built from the same wrong stress.
    This differences the material point on its own first.
    """
    params, state, dstrain = a_yielded_state()
    analytic, yielded = dsigma_dq_closed_form(params, state, dstrain)
    assert yielded, "the state chosen is not plastic, so the derivative is zero"
    q0 = float(state.statev[0])
    sweep = sweep_steps(
        analytic,
        lambda h: (stress_at(params, state, dstrain, q0 + h)
                   - stress_at(params, state, dstrain, q0 - h)) / (2.0 * h),
        steps=STEPS, tolerance=1e-8, scale=max(abs(q0), 1e-6))
    assert not sweep.flat, sweep.verdict()
    assert sweep.converged, (
        f"the closed-form dsigma/dq does not differentiate the stress update: "
        f"{sweep.verdict()}")


def test_the_assembled_state_derivative_is_the_derivative_of_the_assembly():
    """dR/dq on one element against a centred difference of R itself.

    The state is perturbed at ONE integration point, because that is what a
    state is: perturbing all eight at once would pass with the block structure
    thrown away.
    """
    params, state, dstrain = a_yielded_state()
    analytic, _yielded = dsigma_dq_closed_form(params, state, dstrain)
    q0 = float(state.statev[0])

    dsigma_dq = np.zeros((8, 6, 1))
    for ip in range(8):
        dsigma_dq[ip, :, 0] = analytic
    assembled = element_dR_dq(UNIT_CUBE, dsigma_dq)
    assert assembled.shape == (24, 8)

    def residual_with(eqplas_at_ip):
        stresses = np.array([stress_at(params, state, dstrain, q)
                             for q in eqplas_at_ip])
        return element_internal_force_small_strain(UNIT_CUBE, stresses)

    for ip in (0, 3, 7):
        column = assembled[:, state_column(ip, 0, 1)]

        def difference(h, ip=ip):
            plus, minus = np.full(8, q0), np.full(8, q0)
            plus[ip] += h
            minus[ip] -= h
            return (residual_with(plus) - residual_with(minus)) / (2.0 * h)

        sweep = sweep_steps(column, difference, steps=STEPS, tolerance=AGREES,
                            scale=max(abs(q0), 1e-6))
        assert not sweep.flat, f"point {ip}: {sweep.verdict()}"
        assert sweep.converged, (
            f"dR/dq at integration point {ip} is not the derivative of the "
            f"assembled residual: {sweep.verdict()}\n"
            f"worst components: {sweep.worst_components()}")


def test_a_state_at_one_point_moves_that_point_and_no_other():
    """The block structure, asserted rather than assumed.

    ``dR/dq`` built as though the state were shared would have every column
    equal to the sum over points. It is not: the column of point ``k`` is the
    integral over point ``k`` alone.
    """
    params, state, dstrain = a_yielded_state()
    analytic, _ = dsigma_dq_closed_form(params, state, dstrain)
    dsigma_dq = np.zeros((8, 6, 1))
    dsigma_dq[3, :, 0] = analytic                    # one point only
    assembled = element_dR_dq(UNIT_CUBE, dsigma_dq)

    moved = [ip for ip in range(8)
             if np.max(np.abs(assembled[:, state_column(ip, 0, 1)])) > 0.0]
    assert moved == [3], (
        f"a state at point 3 moved the columns of points {moved}; a state "
        f"variable is not shared between integration points")


def test_a_state_derivative_shaped_like_a_parameter_is_refused():
    """The shape that would silently make a state into a parameter."""
    with pytest.raises(StateSensitivityError) as raised:
        element_dR_dq(UNIT_CUBE, np.zeros((6, 1)))
    assert "is not a state derivative" in str(raised.value)


# --------------------------------------------------------------------------- #
# a mesh
# --------------------------------------------------------------------------- #
def a_two_element_mesh():
    coords = []
    for k in (0.0, 1.0):
        for j in (0.0, 1.0):
            for i in (0.0, 1.0, 2.0):
                coords.append((i, j, k))
    node_ids = list(range(1, len(coords) + 1))
    index = {(i, j, k): 1 + i + 3 * j + 6 * k
             for k in range(2) for j in range(2) for i in range(3)}
    connectivity = []
    for e in range(2):
        connectivity.append((e + 1, [
            index[(e, 0, 0)], index[(e + 1, 0, 0)], index[(e + 1, 1, 0)],
            index[(e, 1, 0)], index[(e, 0, 1)], index[(e + 1, 0, 1)],
            index[(e + 1, 1, 1)], index[(e, 1, 1)]]))
    return node_ids, np.array(coords, dtype=float), connectivity


def test_the_global_state_derivative_lands_only_on_the_element_that_owns_it():
    """Scattered across two elements: a state moves its own element's nodes.

    The two elements share a face, so four of the twelve nodes belong to both
    and eight belong to one. A state at an element's integration point may
    move the shared nodes -- they carry that element's force -- and may not
    move the four nodes of the other element that it does not touch.
    """
    params, state, dstrain = a_yielded_state()
    analytic, _ = dsigma_dq_closed_form(params, state, dstrain)
    node_ids, coords, connectivity = a_two_element_mesh()

    per_element = {}
    for eid, _conn in connectivity:
        block = np.zeros((8, 6, 1))
        if eid == 1:
            block[5, :, 0] = analytic                # one point of element 1
        per_element[eid] = block

    dR_dq, node_index, columns = assemble_dR_dq(node_ids, coords, connectivity,
                                                per_element)
    assert dR_dq.shape == (3 * len(node_ids), 2 * 8 * 1)
    column = columns.column_of(1, 5, 0)
    assert columns.describe(column) == {"element": 1, "ip": 5, "state": 0}

    touched = {nid for eid, conn in connectivity if eid == 1 for nid in conn}
    moved_rows = np.where(np.abs(dR_dq[:, column]) > 0.0)[0]
    moved_nodes = {node for node, i in node_index.items()
                   if any(3 * i <= row < 3 * i + 3 for row in moved_rows)}
    assert moved_nodes <= touched, (
        f"a state at element 1 moved nodes {sorted(moved_nodes - touched)}, "
        f"which element 1 does not touch")
    assert moved_nodes, "the state moved nothing at all"

    other = np.zeros(dR_dq.shape[0], dtype=bool)
    for nid in set(node_ids) - touched:
        i = node_index[nid]
        other[3 * i:3 * i + 3] = True
    assert np.all(dR_dq[other, :][:, column] == 0.0)


# --------------------------------------------------------------------------- #
# and the layer it settles
# --------------------------------------------------------------------------- #
def test_a_state_probe_turns_the_state_layer_from_unknown_into_a_verdict():
    """The attribution's third answer, and what removes it.

    Without a dR/dq the state layer reports NOT ESTABLISHED and says what
    would settle it. With one, it reports a plateau. The two must not read the
    same, because "nothing was checked" and "checked and clean" are different
    facts.
    """
    fixture = hex_fixtures()[0]
    bare = diagnose(fixture,
                    assemble=lambda s: element_internal_force_small_strain(
                        UNIT_CUBE, s),
                    reference=lambda s: element_internal_force_small_strain(
                        UNIT_CUBE, s))
    assert "state_sensitivity" in bare.not_established
    assert bare.finding("state_sensitivity").would_establish
    assert not bare.complete and bare.ok

    params, state, dstrain = a_yielded_state()
    analytic, _ = dsigma_dq_closed_form(params, state, dstrain)
    q0 = float(state.statev[0])
    dsigma_dq = np.zeros((8, 6, 1))
    dsigma_dq[2, :, 0] = analytic
    column = element_dR_dq(UNIT_CUBE, dsigma_dq)[:, state_column(2, 0, 1)]

    def difference(h):
        def force(shift):
            stresses = np.array([
                stress_at(params, state, dstrain,
                          q0 + (shift if ip == 2 else 0.0))
                for ip in range(8)])
            return element_internal_force_small_strain(UNIT_CUBE, stresses)
        return (force(h) - force(-h)) / (2.0 * h)

    probe = StateSensitivityProbe(analytic=column, difference_at=difference,
                                  scale=max(abs(q0), 1e-6),
                                  label="dR/dq at integration point 2",
                                  tolerance=AGREES)
    settled = diagnose(fixture, state_probe=probe,
                       assemble=lambda s: element_internal_force_small_strain(
                           UNIT_CUBE, s),
                       reference=lambda s: element_internal_force_small_strain(
                           UNIT_CUBE, s),
                       steps=STEPS)
    finding = settled.finding("state_sensitivity")
    assert finding.ok, finding.message()
    assert finding.measured["plateau"] >= 3
    assert "state_sensitivity" not in settled.not_established


def test_a_wrong_state_derivative_is_blamed_on_the_state_layer():
    """And the failure says which layer, with the sweep behind it.

    The wrong derivative is wrong by a factor, so its error does not move with
    the step size -- which is the signature the verdict has to name, because a
    flat error is a wrong formula and reads like a converged plateau to
    anything that only asks whether the best number is small.
    """
    fixture = hex_fixtures()[0]
    params, state, dstrain = a_yielded_state()
    analytic, _ = dsigma_dq_closed_form(params, state, dstrain)
    q0 = float(state.statev[0])
    dsigma_dq = np.zeros((8, 6, 1))
    dsigma_dq[2, :, 0] = analytic
    column = element_dR_dq(UNIT_CUBE, dsigma_dq)[:, state_column(2, 0, 1)]

    def difference(h):
        def force(shift):
            stresses = np.array([
                stress_at(params, state, dstrain,
                          q0 + (shift if ip == 2 else 0.0))
                for ip in range(8)])
            return element_internal_force_small_strain(UNIT_CUBE, stresses)
        return (force(h) - force(-h)) / (2.0 * h)

    probe = StateSensitivityProbe(analytic=1.5 * column,
                                  difference_at=difference,
                                  scale=max(abs(q0), 1e-6),
                                  tolerance=AGREES)
    found = diagnose(fixture, state_probe=probe, steps=STEPS)
    assert found.blame == "state_sensitivity", found.report()
    assert found.finding("state_sensitivity").measured["flat"] is True
    assert "wrong formula" in found.report()
