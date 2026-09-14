"""Three answers reach this repository, and only two of them are a verdict.

The assembler's job is to build R from ingredients, two of which come from a
UMAT. If it accepts an ingredient whose verification gate was never measured,
it produces a number that looks exactly like a verified one. So the reading is
three-state here, on the consuming side, and not only where the record was
written: a claim checked only by the tool that makes it is a claim nobody
checked.
"""
from __future__ import annotations

import json

import pytest

from contract_reader import (GATES, SEVENTH, ContractError, gates_of,
                             primal_settled, usable_for_assembly)
from tristate import NOT_ESTABLISHED, Tri, TristateError, all_true, read

STORE_FINGERPRINT = "b0d27ee53c630500"


def _record(**overrides) -> dict:
    base = {
        "contract_version": "1.0.0",
        "identity": {"path": "repo__x/src/umat.f", "sha256": "a" * 64},
        "transform_fingerprint": STORE_FINGERPRINT,
        "terminal": {"state": "fully_verified", "owner": "NONE"},
        "evidence": {name: True for name in GATES},
        "convention": {"voigt_order": ["11", "22", "33", "12", "13", "23"],
                       "shear": "engineering"},
    }
    base.update(overrides)
    return base


def test_the_reader_refuses_to_be_used_as_a_boolean():
    with pytest.raises(TristateError) as exc:
        bool(NOT_ESTABLISHED)
    assert "is_not_established" in str(exc.value)
    with pytest.raises(TristateError):
        bool(Tri(True))


def test_a_missing_gate_and_a_null_gate_are_the_same_answer():
    absent = gates_of({"evidence": {}})
    null = gates_of({"evidence": {name: None for name in GATES}})
    for name in GATES:
        assert absent[name].is_not_established()
        assert null[name].is_not_established()
        assert absent[name] == null[name]


def test_a_record_with_no_evidence_block_is_not_a_verified_record():
    record = _record()
    del record["evidence"]
    answer, why = usable_for_assembly(record)
    assert answer.is_not_established()
    assert not answer.is_true()
    assert "never measured" in why


def test_one_unmeasured_gate_is_enough_to_withhold_the_verdict():
    for missing in GATES:
        evidence = {name: True for name in GATES if name != missing}
        answer, why = usable_for_assembly(_record(evidence=evidence))
        assert not answer.is_true(), missing
        assert answer.is_not_established() or answer.is_false()


def test_a_measured_failure_and_an_unmeasured_gate_are_different_answers():
    failed = usable_for_assembly(_record(
        evidence={**{n: True for n in GATES}, "derivatives_verified": False}))
    unmeasured = usable_for_assembly(_record(
        evidence={n: True for n in GATES if n != "derivatives_verified"}))
    assert failed[0].is_false()
    assert unmeasured[0].is_not_established()
    assert failed[0] != unmeasured[0]
    # A refusal is a finding; an unmeasured entry is a queue item. Merging
    # them puts work that cannot be done beside work nobody has started.
    assert "measured false" in failed[1]
    assert "never measured" in unmeasured[1]


def test_there_are_six_gates_and_the_seventh_is_not_one_of_them():
    assert len(GATES) == 6
    assert SEVENTH not in GATES
    assert GATES == ("abaqus_job_completed", "all_requested_outputs_present",
                     "complete_history_finite", "primal_agreed",
                     "derivatives_verified", "mechanically_informative")


def test_an_unexplained_disagreement_is_never_a_verification():
    """The builds disagreed and nothing says whether a control explained it.
    NOT ESTABLISHED, and NOT ESTABLISHED is not a pass."""
    gates = gates_of({"evidence": {**{n: True for n in GATES},
                                   "primal_agreed": False}})
    settled = primal_settled(gates)
    assert settled.is_not_established()
    assert not settled.is_true()
    answer, _ = usable_for_assembly(_record(
        evidence={**{n: True for n in GATES}, "primal_agreed": False}))
    assert not answer.is_true()


def test_a_disagreement_a_control_explained_is_a_chain_not_a_contradiction():
    """13 of the 55 verified entries in the frozen store have ``primal_agreed``
    measured FALSE, and their verdicts rest on a control having measured why.
    A consumer that read the raw gate alone would refuse all 13."""
    evidence = {**{n: True for n in GATES}, "primal_agreed": False,
                SEVENTH: True}
    gates = gates_of({"evidence": evidence})
    assert primal_settled(gates).is_true()
    answer, _ = usable_for_assembly(_record(evidence=evidence))
    assert answer.is_true()
    # And the raw comparison never moves: it still says the builds differed.
    assert gates["primal_agreed"].is_false()


def test_a_control_that_ran_and_explained_nothing_is_a_refusal():
    evidence = {**{n: True for n in GATES}, "primal_agreed": False,
                SEVENTH: False}
    assert primal_settled(gates_of({"evidence": evidence})).is_false()
    answer, _ = usable_for_assembly(_record(evidence=evidence))
    assert answer.is_false()


def test_a_truthy_string_cannot_reach_a_gate():
    for smuggled in ("true", "false", "no", "0", 1, 0):
        with pytest.raises(TristateError):
            gates_of({"evidence": {"primal_agreed": smuggled}})


def test_an_evidence_block_that_is_not_an_object_is_refused():
    for wrong in (True, "verified", 1, ["ok"]):
        with pytest.raises(ContractError):
            gates_of({"evidence": wrong})


def test_the_conjunction_keeps_the_third_state():
    assert all_true(Tri(True), Tri(False), NOT_ESTABLISHED).is_false()
    assert all_true(Tri(True), NOT_ESTABLISHED).is_not_established()
    assert all_true(Tri(True), Tri(True)).is_true()


def test_an_entry_with_no_tensor_convention_cannot_be_assembled():
    """An assembled residual indexed by a component order nobody stated is a
    number that is right in some components and wrong in others."""
    answer, why = usable_for_assembly(_record(convention={"voigt_order": []}))
    assert answer.is_not_established()
    assert "component order" in why


# ---------------------------------------------------------------------------
# what is actually committed here
# ---------------------------------------------------------------------------
def test_a_frozen_fixture_may_not_carry_a_false_primal_with_no_explanation():
    """Found in committed data, at the current transform generation.

    A fixture can carry all six gates and have ``primal_agreed`` be one of
    them and be FALSE -- the two builds differed. If it does not also carry
    ``primal_difference_explained_by_a_measured_control``, nothing in it says
    whether a control measured the reason, and under this contract that reads
    NOT ESTABLISHED. Not a pass, and not a refusal either.

    The fix is upstream: the exporter writes the six and should carry the
    seventh with them, because a fixture whose own evidence block says a gate
    failed and says nothing about why is a regression baseline a reader cannot
    evaluate. Until it does, this is what stops the assembler resting on one.
    """
    from pathlib import Path

    from contract_paths import FIXTURES
    from contract_reader import primal_settled

    for path in sorted(FIXTURES.glob("*.json")):
        payload = json.loads(path.read_text())
        evidence = (payload.get("finite_history") or {}).get("evidence") or {}
        gates = gates_of({"evidence": evidence})
        if not gates["primal_agreed"].is_false():
            continue
        assert gates[SEVENTH].is_not_established()
        settled = primal_settled(gates)
        assert settled.is_not_established()
        assert not settled.is_true()
        answer, why = usable_for_assembly({
            "terminal": {"state": "fully_verified", "owner": "NONE"},
            "evidence": evidence,
            "convention": {"voigt_order": ["11", "22", "33", "12", "13", "23"]}})
        assert not answer.is_true(), path.name


def test_no_committed_fixture_hides_a_gate_behind_a_truthy_value():
    """Every gate in every committed fixture is true, false or null. A string
    or a number would be truthy, and would pass a reader that tested it."""
    from contract_paths import FIXTURES

    for path in sorted(FIXTURES.glob("*.json")):
        payload = json.loads(path.read_text())
        evidence = (payload.get("finite_history") or {}).get("evidence") or {}
        for name, value in evidence.items():
            assert value is None or isinstance(value, bool), \
                f"{path.name}: {name} is {value!r}"
        gates_of({"evidence": evidence})       # raises on anything else


def test_the_reason_names_the_primal_chain_and_not_just_the_raw_gate():
    """When the only unresolved thing is an unexplained disagreement, the
    reason has to say that -- "a gate was never measured:" with nothing after
    it tells a reader nothing and sends them to the wrong place."""
    evidence = {**{n: True for n in GATES}, "primal_agreed": False}
    answer, why = usable_for_assembly(_record(evidence=evidence))
    assert answer.is_not_established()
    assert "primal_agreed" in why
    assert "unexplained disagreement is not a verification" in why
    assert not why.rstrip().endswith(":")

    refused = {**evidence, SEVENTH: False}
    answer, why = usable_for_assembly(_record(evidence=refused))
    assert answer.is_false()
    assert "primal_agreed" in why and "did not explain" in why
