"""The diagnostic layers, broken on purpose, and the verdict that follows.

An attribution nobody has watched fail is a list of strings. Every layer below
is given a fault that belongs to it and nothing else, and the diagnosis has to
name that layer -- not the one before it, not the one after it. The faults are
the real failure modes: a NaN, two builds that disagree, a tangent that is not
a derivative, a state derivative wrong by a factor, a shear convention, a
quadrature weight, an assembly that does not match an independent one, a
partition that leaks, a DOF map with a collision.

And the answers that are not verdicts. A layer this diagnosis was given
nothing to check reports NOT ESTABLISHED and says what would settle it, which
is a different fact from holding, and the tests here hold those two apart --
because the whole value of the attribution is that a reader can tell "checked
and clean" from "nobody looked".
"""
from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from residual_core.core.deck_replay import replay_deck  # noqa: E402
from residual_core.core.where_it_went_wrong import (  # noqa: E402
    FAILS, HOLDS, LAYERS, NOT_ESTABLISHED, OWNER, diagnose)
from residual_core.formulations.c3d8_kernel import (  # noqa: E402
    element_internal_force_small_strain)
from residual_core.materials.verified_fixture import (  # noqa: E402
    CURRENT_TRANSFORM_FINGERPRINT, EVIDENCE_GATES, FixtureError, load,
    tangent_convention)
from verified_fixtures import (FIXTURES, UNIT_CUBE, all_fixtures,  # noqa: E402
                               fixture_named, hex_fixtures)


def assembled(stress):
    return element_internal_force_small_strain(UNIT_CUBE, stress)


def run(fixture, **kwargs):
    kwargs.setdefault("assemble", assembled)
    kwargs.setdefault("reference", assembled)
    return diagnose(fixture, **kwargs)


# --------------------------------------------------------------------------- #
# the shape of the answer
# --------------------------------------------------------------------------- #
def test_the_layers_include_parameter_sensitivity_and_element_jacobian():
    assert LAYERS == ("umat", "transformation", "constitutive_derivative",
                      "state_sensitivity", "parameter_sensitivity", "mapping",
                      "element_integration", "residual_assembly",
                      "element_jacobian", "boundary_conditions",
                      "global_dof_mapping")
    assert set(OWNER) == set(LAYERS)
    for layer, owner in OWNER.items():
        assert owner and owner[0].islower(), layer


def test_every_verdict_carries_the_numbers_it_was_reached_on():
    """A layer named without evidence is worse than no attribution."""
    found = run(hex_fixtures()[0])
    for finding in found.findings:
        if finding.status == NOT_ESTABLISHED:
            assert finding.would_establish, finding.stage
            continue
        assert finding.detail, finding.stage
        assert finding.measured, (
            f"{finding.stage} reports {finding.status} with no measurement "
            f"behind it")
        assert any(isinstance(v, (int, float)) for v in finding.measured.values()), (
            f"{finding.stage} carries no number: {finding.measured}")


def test_not_established_is_not_a_pass_and_says_so_in_the_report():
    found = run(hex_fixtures()[0])
    assert "state_sensitivity" in found.not_established
    assert found.ok, "an unestablished layer is not a failure"
    assert not found.complete, "and it is not a clean bill either"
    assert "NOT ESTABLISHED" in found.report()
    assert "NOT established: state_sensitivity" in found.report()
    finding = found.finding("state_sensitivity")
    assert finding.status == NOT_ESTABLISHED and not finding.ok


def test_the_whole_diagnosis_serialises_for_a_caller_that_is_not_a_person():
    """The GUI reads this, not the prose."""
    payload = run(hex_fixtures()[0]).as_dict()
    assert json.loads(json.dumps(payload, default=float))
    assert set(payload) == {"ok", "complete", "blame", "not_established",
                            "findings"}
    assert [f["layer"] for f in payload["findings"]] == list(LAYERS)
    for finding in payload["findings"]:
        assert finding["status"] in (HOLDS, FAILS, NOT_ESTABLISHED)


# --------------------------------------------------------------------------- #
# provenance: the layer before the layers
# --------------------------------------------------------------------------- #
def test_a_fixture_from_an_older_transformation_is_refused_by_the_loader(tmp_path):
    """The refusal that would have saved this repository a session of work.

    Every fixture committed before 2026-09-14 carried ff94800b1884bcc0 and
    nothing read it, so every comparison was against a transformation two
    builds old and each one of them passed.
    """
    payload = json.loads((FIXTURES / "isotropic-elasticity--f7eb90376a.json")
                         .read_text(encoding="utf-8"))
    payload["transform_fingerprint"] = "ff94800b1884bcc0"
    stale = tmp_path / "stale.json"
    stale.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(FixtureError) as raised:
        load(stale)
    said = str(raised.value)
    assert "ff94800b1884bcc0" in said and CURRENT_TRANSFORM_FINGERPRINT in said
    assert "somebody else's run" in said
    # and reading it deliberately, with a reason, is still possible
    assert load(stale, fingerprint=None).transform_fingerprint == \
        "ff94800b1884bcc0"


def test_a_fixture_with_no_fingerprint_at_all_is_refused(tmp_path):
    payload = json.loads((FIXTURES / "umat--b8fa38353e.json")
                         .read_text(encoding="utf-8"))
    payload.pop("transform_fingerprint")
    bare = tmp_path / "bare.json"
    bare.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(FixtureError) as raised:
        load(bare)
    assert "none recorded" in str(raised.value)


def test_every_committed_fixture_carries_the_gates_the_pipeline_measured():
    for fixture in all_fixtures():
        assert set(fixture.evidence) >= set(EVIDENCE_GATES), fixture.path.name
        assert set(fixture.gates_true) | set(fixture.gates_not_true) == \
            set(EVIDENCE_GATES)


# --------------------------------------------------------------------------- #
# which reading of DDSDDE -- measured per material, not per flag
# --------------------------------------------------------------------------- #
def test_the_tangent_reading_is_decided_by_the_numbers_and_not_by_nlgeom():
    """Two fixtures, the same NLGEOM=YES, opposite readings.

    irfancn's elastic UMAT is a small-strain routine handed a logarithmic
    strain: its DDSDDE satisfies the plain reading. AlexanderJFDR's
    neo-Hookean under the same flag satisfies the Jaumann/Kirchhoff one. An
    assembler that chose by the flag would build the wrong stiffness for one
    of them, and the error is 8.7e-03 relative -- big enough to matter and
    small enough to be argued about.
    """
    plain = tangent_convention(fixture_named("umat_elastic--7e9bb4c291.json"))
    jaumann = tangent_convention(fixture_named("neohookean_umat--10759f1ffd.json"))
    assert plain.reading == "material" and jaumann.reading == "jaumann"
    assert plain.material_error < 1e-14 < 1e-3 < plain.jaumann_error
    assert jaumann.jaumann_error < 1e-6 < 1e-3 < jaumann.material_error
    assert plain.separation > 1e6 and jaumann.separation > 1e4
    for convention in (plain, jaumann):
        assert convention.volumetric > 0.0, (
            "an isochoric window cannot separate the readings and must not "
            "claim to")


def test_an_isochoric_window_refuses_to_choose_a_reading():
    """The two readings differ by sigma tr(Deps). With a zero trace they are
    the same formula, and saying which one holds would be inventing a fact."""
    fixture = fixture_named("neohookean_umat--10759f1ffd.json")
    flattened = dataclasses.replace(fixture)
    flattened.converted = [
        dataclasses.replace(record, dstrain=np.array(
            [0.0, 0.0, 0.0, *record.dstrain[3:]]))
        for record in fixture.converted]
    convention = tangent_convention(flattened)
    assert convention.reading is None
    assert "isochoric" in convention.why


def test_the_pipelines_own_gate_and_this_repositorys_layer_agree(capsys):
    """The viscoelastic case the store recorded primal_agreed=False for.

    Two independent readings of the same run: the pipeline's gate, recorded
    in the fixture, and this repository's transformation layer, computed from
    the carried numbers. They must reach the same verdict, and the layer must
    put a number on it.
    """
    fixture = fixture_named("umat_viscoelastic--c6ae96a734.json")
    assert fixture.evidence["primal_agreed"] is False
    assert fixture.gates_not_true == ("primal_agreed",), (
        "this case is the one that failed exactly one gate; if that changed, "
        "the test below is measuring something else")
    found = run(fixture)
    assert found.blame == "transformation", found.report()
    measured = found.finding("transformation").measured
    assert 1e-9 < measured["stress"] < 1e-6, (
        f"the pipeline recorded a disagreement and this layer measured "
        f"{measured['stress']:.3e}; the two readings have to be the same size")


# --------------------------------------------------------------------------- #
# one fault per layer
# --------------------------------------------------------------------------- #
def test_a_shear_convention_in_the_tangent_is_blamed_on_the_mapping_layer():
    """Halving the shear block of the tangent is exactly the tensorial/
    engineering mistake, and it is a mapping fault rather than a wrong
    tangent."""
    fixture = fixture_named("isotropic-elasticity--f7eb90376a.json")
    broken = dataclasses.replace(fixture)
    broken.converted = []
    for record in fixture.converted:
        tangent = record.tangent.copy()
        tangent[3:, 3:] *= -1.0            # no stable material has this
        broken.converted.append(dataclasses.replace(record, tangent=tangent))
    found = run(broken)
    assert found.blame == "mapping", found.report()
    assert "shear diagonal entry is not positive" in found.report()


def test_a_window_with_no_shear_leaves_the_shear_convention_unestablished():
    """The honest answer where the numbers are silent.

    Every committed fixture's window sits in its deck's uniaxial step, so a
    factor of two on the shear components changes nothing in it. The mapping
    layer says so instead of reporting agreement.
    """
    found = run(hex_fixtures()[0])
    mapping = found.finding("mapping")
    assert mapping.status == NOT_ESTABLISHED, mapping.message()
    assert mapping.measured["shear_exercised"] < 1e-12
    assert "factor of two" in mapping.detail
    assert "shear step" in mapping.would_establish


def test_a_quadrature_that_does_not_integrate_is_blamed_on_the_element_layer():
    """A weight that is 2% wrong. The divergence theorem is not a tolerance."""
    from residual_core.formulations.c3d8_kernel import _GaussRule, ABAQUS_C3D8_GAUSS

    wrong = _GaussRule(ABAQUS_C3D8_GAUSS.points,
                       ABAQUS_C3D8_GAUSS.weights * 1.02)
    found = run(hex_fixtures()[0], gauss=wrong)
    assert found.blame == "element_integration", found.report()
    assert found.finding("element_integration").measured[
        "uniform_stress_vs_tractions"] > 1e-3


def test_an_element_geometry_that_is_degenerate_is_blamed_on_the_element_layer():
    squashed = UNIT_CUBE.copy()
    squashed[4:, 2] = 0.0                  # zero thickness
    found = run(hex_fixtures()[0], coordinates=squashed)
    assert found.blame == "element_integration", found.report()


def test_a_partition_that_leaks_is_blamed_on_the_boundary_layer():
    fixture = hex_fixtures()[0]
    replay = replay_deck(fixture.deck, step=1)
    replay.model.boundaries = []           # a step that holds nothing
    replay.prescribed = {}
    found = run(fixture, replay=replay)
    assert found.blame == "boundary_conditions", found.report()
    assert "null space" in found.report()


def test_a_dof_map_with_a_collision_is_blamed_on_the_mapping_of_dofs():
    fixture = hex_fixtures()[0]
    replay = replay_deck(fixture.deck, step=1)

    class Collides:
        """A map that sends two different unknowns to one row."""

        def __init__(self, inner):
            self._inner = inner
            self.node_ids = inner.node_ids
            self.ndof = inner.ndof

        def __getattr__(self, name):
            return getattr(self._inner, name)

        def dof_index(self, nid, local):
            return 0 if (nid, local) == (8, 3) else \
                self._inner.dof_index(nid, local)

    replay.dof_manager = Collides(replay.dof_manager)
    found = run(fixture, replay=replay)
    assert found.blame == "global_dof_mapping", found.report()
    assert "cannot share a row" in found.report()


def test_without_a_deck_the_last_two_layers_are_unknown_rather_than_clean():
    """A fixture with no deck leaves the boundary and DOF layers unchecked,
    and the report says which."""
    fixture = dataclasses.replace(hex_fixtures()[0], deck="")
    found = run(fixture, use_deck=True)
    assert set(found.not_established) >= {"boundary_conditions",
                                          "global_dof_mapping"}
    assert found.ok and not found.complete
    for layer in ("boundary_conditions", "global_dof_mapping"):
        assert found.finding(layer).would_establish


def test_the_run_stops_at_the_first_layer_that_fails():
    """Everything downstream of a NaN is meaningless, and reporting nine
    verdicts would bury the one that matters."""
    fixture = hex_fixtures()[0]
    broken = dataclasses.replace(fixture)
    broken.original = [dataclasses.replace(r) for r in fixture.original]
    broken.original[-1] = dataclasses.replace(
        broken.original[-1],
        stress=np.array([float("nan"), *fixture.original[-1].stress[1:]]))
    found = run(broken)
    assert found.blame == "umat"
    assert len(found.findings) == 1
    assert OWNER["umat"] in found.report()
