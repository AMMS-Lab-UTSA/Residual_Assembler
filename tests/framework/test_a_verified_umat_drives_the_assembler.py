"""The assembler consumes what the pipeline actually produced, not a stand-in.

Its tests used to supply stress and tangent from a reference model written for
the purpose. That tests the assembly and nothing about the bridge: a fixture
written to make the assembler pass cannot also be evidence that the assembler
reads what the UMAT pipeline emits.

These read a corpus case that VERIFIED in Abaqus -- the original UMAT ran, the
OTI-converted build ran on the same deck, their stress and state histories
agreed over the whole path, and the converted build's tangent agreed with a
finite difference of the original at several states. What the fixture carries
is the numbers and the identity of the file they came from; the file itself
stays out, because most of the corpus is not redistributable.

Two things are being checked. That the residual and its derivative assembled
from those numbers are right -- against finite differences of the assembly
itself, and against an independently written integration. And that when
something is wrong, the failure says WHICH of the five stages it came from,
because "the residual is out by 12%" tells nobody which of five people to ask.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from residual_core.core.where_it_went_wrong import (OWNER, STAGES,  # noqa: E402
                                                    diagnose)
from residual_core.formulations.c3d8_kernel import (  # noqa: E402
    ABAQUS_C3D8_GAUSS, b_matrix_reference, element_internal_force_small_strain,
    element_tangent)
from residual_core.materials.verified_fixture import (  # noqa: E402
    FixtureError, SCHEMA, check_conventions, load)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verified_fixtures import (FIXTURES, UNIT_CUBE,  # noqa: E402
                               all_fixtures, hex_fixtures)


def fixtures():
    """Only the hexahedra: the checks below integrate ``B^T sigma`` on a unit
    cube, and a three-component plane-stress tensor is not a smaller case of
    that -- it is a different element."""
    return hex_fixtures()


def independent_internal_force(Xe, stress_at_ip):
    """B^T sigma integrated here, written out rather than called.

    A second implementation, so that stage five is checked against something
    other than itself. Deliberately the long way round.
    """
    out = np.zeros(24)
    for index, (point, weight) in enumerate(zip(ABAQUS_C3D8_GAUSS.points,
                                                ABAQUS_C3D8_GAUSS.weights)):
        B, detJ = b_matrix_reference(Xe, point)
        out = out + (B.T @ np.asarray(stress_at_ip[index], dtype=float)) * detJ * weight
    return out


def test_every_fixture_declares_what_it_is():
    for fixture in fixtures():
        assert fixture.source_id, fixture.path
        assert fixture.material_provenance, (
            f"{fixture.path.name} carries constants with no stated origin, "
            f"which reads as established when nothing established it")
        assert fixture.verification.get("states_checked"), fixture.path.name
        assert fixture.deck, "the experiment that produced these numbers"


def test_a_fixture_of_another_shape_is_refused(tmp_path: Path):
    """A fixture whose shape is not the shape this reads would be interpreted
    wrongly rather than refused."""
    bad = tmp_path / "bad.json"
    bad.write_text('{"schema": "something-else/9"}', encoding="utf-8")
    with pytest.raises(FixtureError) as raised:
        load(bad)
    assert SCHEMA in str(raised.value)


def test_a_fixture_that_contradicts_its_own_tensor_size_is_refused(tmp_path: Path):
    import json
    source = fixtures()[0]
    payload = json.loads(source.path.read_text(encoding="utf-8"))
    payload["material_point"]["ntens"] = 4
    bad = tmp_path / "wrong_ntens.json"
    bad.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(FixtureError) as raised:
        load(bad)
    assert "NTENS=4" in str(raised.value)


def test_the_conventions_the_fixture_states_hold_in_its_numbers():
    for fixture in fixtures():
        assert check_conventions(fixture) == [], fixture.path.name


def test_the_internal_force_matches_an_independent_integration():
    for fixture in fixtures():
        stress = np.tile(fixture.original[-1].stress, (8, 1))
        mine = element_internal_force_small_strain(UNIT_CUBE, stress)
        theirs = independent_internal_force(UNIT_CUBE, stress)
        assert np.allclose(mine, theirs, rtol=0, atol=1e-12), fixture.path.name


def test_a_uniform_stress_leaves_no_net_force():
    """Equilibrium: a constant stress field on a free element produces
    internal forces that sum to zero over the nodes.

    Judged relative to the forces being summed, not against a fixed absolute
    number. The corpus fixtures span eight orders of magnitude of stress --
    irfancn__Abaqus-UMAT-elastic reaches 1.3e+08 where the bundled J2 control
    is at 1.7e+03 -- and a cancellation of 1.3e+08 against itself leaves
    7.5e-09 of float64 round-off, which an absolute 1e-9 reads as a violated
    equilibrium and is in fact 5.7e-17 of the quantity that cancelled.
    """
    for fixture in fixtures():
        stress = np.tile(fixture.original[-1].stress, (8, 1))
        force = element_internal_force_small_strain(UNIT_CUBE, stress).reshape(8, 3)
        scale = max(float(np.max(np.abs(force))), 1.0)
        assert np.max(np.abs(force.sum(axis=0))) <= 1e-12 * scale, (
            f"{fixture.path.name}: net force "
            f"{force.sum(axis=0)} against forces of {scale:.3e}")


def test_the_element_tangent_is_the_derivative_of_the_internal_force():
    """The assembled dR/du against a centred difference of the assembly.

    The material tangent comes from the verified fixture, so what is being
    checked is this repository's integration of it -- not the UMAT's, which
    Abaqus and the finite-difference replay already checked.
    """
    for fixture in fixtures():
        record = fixture.converted[-1]
        if record.tangent is None or fixture.ntens != 6:
            continue
        D = np.tile(record.tangent, (8, 1, 1))
        K = element_tangent(UNIT_CUBE, np.zeros(24), D, mode="small")

        def force(displacement):
            strain_at_ip = []
            for point in ABAQUS_C3D8_GAUSS.points:
                B, _detJ = b_matrix_reference(UNIT_CUBE, point)
                strain_at_ip.append(B @ displacement)
            stress = np.array([record.tangent @ e for e in strain_at_ip])
            return element_internal_force_small_strain(UNIT_CUBE, stress)

        step = 1e-7
        numerical = np.zeros((24, 24))
        for dof in range(24):
            plus, minus = np.zeros(24), np.zeros(24)
            plus[dof], minus[dof] = step, -step
            numerical[:, dof] = (force(plus) - force(minus)) / (2.0 * step)
        scale = float(np.max(np.abs(K)))
        assert np.max(np.abs(K - numerical)) / scale < 1e-6, (
            f"{fixture.path.name}: the assembled tangent is not the derivative "
            f"of the assembled force")


def test_a_rigid_translation_produces_no_force():
    for fixture in fixtures():
        record = fixture.converted[-1]
        if record.tangent is None or fixture.ntens != 6:
            continue
        translation = np.tile(np.array([0.3, -0.2, 0.1]), 8)
        strain = []
        for point in ABAQUS_C3D8_GAUSS.points:
            B, _detJ = b_matrix_reference(UNIT_CUBE, point)
            strain.append(B @ translation)
        assert np.allclose(strain, 0.0, atol=1e-12), (
            "a rigid translation is not a strain")


# ---------------------------------------------------------------------------
# and when something is wrong, which of the five stages was it?
# ---------------------------------------------------------------------------
def _diagnose(fixture):
    return diagnose(
        fixture,
        assemble=lambda stress: element_internal_force_small_strain(
            UNIT_CUBE, stress),
        reference=lambda stress: independent_internal_force(UNIT_CUBE, stress))


#: The one committed fixture whose material is not linear over an increment,
#: and what the diagnosis finds in it. Recorded rather than excluded.
#:
#: ``constitutive_derivative`` asks whether the reported tangent predicts the
#: reported stress increment -- a SECANT across a finite increment. For an
#: elastic or hyperelastic material driven in small steps the secant and the
#: tangent agree, which is why the other nine pass it. J2 is the first
#: committed fixture with a real nonlinearity, and its worst disagreement,
#: 3.102e-01, is at INCREMENT 2: the increment in which it yields, where
#: EQPLAS goes from exactly 0 to 2.478519e-04. The tangent at the start of
#: that increment is the elastic one and the stress change across it is
#: elastoplastic, so no tangent evaluated at either end predicts it.
#:
#: That is a limitation of the check, not a defect in the material: this
#: project's standing rule is that a derivative is never evaluated across
#: yielding, damage initiation or any other nonsmooth transition, and this
#: check has no notion of one. The fixture reaches `transformation: holds` at
#: 3.210e-16 in stress and 2.158e-16 in state over all 35 increments before it
#: stops, so what is established about it is established.
#:
#: Written down as a measured fact so that it cannot be mistaken for a pass
#: and cannot be quietly lost: if the check learns about nonsmooth increments,
#: or if the numbers move, this fails and someone reads it.
NONLINEAR_OVER_AN_INCREMENT = {
    "bundled__generic_ps/src/j2_props.f": "constitutive_derivative",
}


def test_a_verified_fixture_passes_every_stage():
    for fixture in fixtures():
        found = _diagnose(fixture)
        expected = NONLINEAR_OVER_AN_INCREMENT.get(fixture.source_id)
        if expected is None:
            assert found.ok, found.report()
        else:
            assert not found.ok, (
                f"{fixture.source_id} now passes every stage. If the "
                f"constitutive_derivative check learned to skip the increment "
                f"a material yields in, delete its entry from "
                f"NONLINEAR_OVER_AN_INCREMENT and say so.")
            assert found.blame == expected, found.report()
            # and everything before the stage it stops at really did hold
            for finding in found.findings:
                if finding.stage == expected:
                    break
                assert finding.status == "holds", found.report()
            # The diagnosis stops at the first failure, so it reports the
            # stages up to and including that one and no more -- reporting
            # later stages it never reached would be claiming they were
            # checked.
            reached = [f.stage for f in found.findings]
            assert reached == list(STAGES)[:len(reached)]
            assert reached[-1] == expected
            continue
        assert [f.stage for f in found.findings] == list(STAGES)


def test_the_yield_increment_is_what_the_derivative_check_trips_on():
    """Named, so the reason J2 stops where it does is not guessed at later.

    The worst disagreement is at the increment where EQPLAS leaves zero. A
    secant taken across a yield point is not a derivative of anything.
    """
    j2 = [f for f in fixtures()
          if f.source_id in NONLINEAR_OVER_AN_INCREMENT]
    if not j2:
        import pytest as _pytest
        _pytest.skip("the J2 fixture is not committed here")
    fixture = j2[0]
    found = _diagnose(fixture)
    measured = next(f.measured for f in found.findings
                    if f.stage == "constitutive_derivative")
    worst_at = int(measured["increment"])

    plastic = [r.state[0] for r in fixture.original]
    assert plastic[0] == 0.0, "the first increment is elastic"
    first_yield = next(i for i, q in enumerate(plastic, start=1) if q > 0.0)
    assert worst_at == first_yield, (
        f"the derivative check's worst increment is {worst_at} and the "
        f"material first yields at {first_yield}; if those have come apart, "
        f"the explanation recorded in NONLINEAR_OVER_AN_INCREMENT no longer "
        f"describes what is happening")


def test_a_stress_that_is_not_a_number_is_blamed_on_the_umat():
    fixture = fixtures()[0]
    fixture.original[-1].stress[0] = float("nan")
    found = _diagnose(fixture)
    assert found.blame == "umat"
    assert "not numbers" in found.report()
    assert OWNER["umat"] in found.report()


def test_a_converted_build_that_disagrees_is_blamed_on_the_transformation():
    fixture = fixtures()[0]
    fixture.converted[-1].stress[0] *= 1.01
    found = _diagnose(fixture)
    assert found.blame == "transformation"
    assert "not the same model" in found.report()


def test_a_tangent_that_does_not_differentiate_its_stress_is_blamed_on_it():
    import dataclasses

    fixture = fixtures()[0]
    fixture.converted = [
        dataclasses.replace(record, tangent=record.tangent * 2.0)
        if record.tangent is not None else record
        for record in fixture.converted]
    found = _diagnose(fixture)
    assert found.blame == "constitutive_derivative"
    assert "is not a tangent" in found.report()


def test_an_assembly_that_disagrees_with_an_independent_one_is_blamed_on_it():
    fixture = fixtures()[0]
    found = diagnose(
        fixture,
        assemble=lambda stress: element_internal_force_small_strain(
            UNIT_CUBE, stress) * 1.5,
        reference=lambda stress: independent_internal_force(UNIT_CUBE, stress))
    assert found.blame == "residual_assembly"
    assert "this is the integration" in found.report()


def test_the_blame_stops_at_the_first_stage_that_fails():
    """Everything downstream of a NaN is meaningless, and reporting all five
    would bury the one that matters."""
    fixture = fixtures()[0]
    fixture.original[-1].stress[0] = float("nan")
    fixture.converted[-1].stress[0] *= 1.5
    found = _diagnose(fixture)
    assert found.blame == "umat"
    assert len(found.findings) == 1


def test_every_stage_names_who_fixes_it():
    assert set(OWNER) == set(STAGES)
    for stage, owner in OWNER.items():
        assert owner and owner[0].islower(), stage
