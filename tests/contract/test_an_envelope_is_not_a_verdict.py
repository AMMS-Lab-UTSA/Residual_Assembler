"""237 calls succeeded. 42 materials are true on all six gates.

``ServiceResult.ok`` says the producing service RAN. It is true of every entry
the service managed to look at, including every entry it looked at and
refused. Reading it as the answer offered all 237 corpus materials to this
repository as verified.

A boolean that is true 237 times out of 237 is not a filter; it is the absence
of one. So the success of the call and the verdict about the subject carry
different names AND different types here: the second is three-state and cannot
be used as a condition at all, because two booleans side by side are two
things that look alike and a reader reaches for whichever looks like success.
"""
from __future__ import annotations

import json

import pytest

from contract_paths import SCHEMAS
from contract_reader import (FORBIDDEN_AS_A_VERDICT, GATES, SUCCESS_FIELD,
                             VERDICT_FIELD, ContractError, usable_for_assembly,
                             verdict_of)
from tristate import NOT_ESTABLISHED, Tri, TristateError


def test_the_two_facts_have_different_names():
    assert SUCCESS_FIELD == "call_succeeded"
    assert VERDICT_FIELD == "verdict"
    assert SUCCESS_FIELD != "ok"


def test_ok_is_never_read_as_a_verdict():
    for shape in ({"ok": True}, {"ok": True, "outcome": "refused"},
                  {"success": True}, {"passed": True}, {"status": "done"}):
        with pytest.raises(ContractError) as exc:
            verdict_of(shape)
        assert VERDICT_FIELD in str(exc.value)
    assert "ok" in FORBIDDEN_AS_A_VERDICT
    assert "RAN" in FORBIDDEN_AS_A_VERDICT["ok"]


def test_the_refusal_names_the_damage_it_prevents():
    with pytest.raises(ContractError) as exc:
        verdict_of({"ok": True})
    assert "237" in str(exc.value) and "42" in str(exc.value)


def test_a_payload_with_no_verdict_raises_rather_than_shrugging():
    """Answering NOT ESTABLISHED would hide the misreading: the caller asked a
    question this payload does not answer, and needs telling."""
    with pytest.raises(ContractError) as exc:
        verdict_of({"service": "verify"})
    assert "nothing that could be mistaken for one" in str(exc.value)


def test_a_stated_verdict_is_read_three_state():
    assert verdict_of({VERDICT_FIELD: True}).is_true()
    assert verdict_of({VERDICT_FIELD: False}).is_false()
    assert verdict_of({VERDICT_FIELD: None}).is_not_established()
    with pytest.raises(ContractError):
        verdict_of({VERDICT_FIELD: "true"})


def test_a_verdict_cannot_be_used_as_a_condition():
    with pytest.raises(TristateError):
        bool(verdict_of({VERDICT_FIELD: True}))
    with pytest.raises(TristateError):
        bool(NOT_ESTABLISHED)


def test_the_schema_refuses_an_envelope_carrying_ok():
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((SCHEMAS / "umat_contract_v1.schema.json").read_text())
    validator = jsonschema.Draft7Validator(schema)

    def record(envelope):
        return {"contract_version": "2.0.0",
                "identity": {"path": "repo__x/u.f", "sha256": "a" * 64},
                "transform_fingerprint": "b0d27ee53c630500",
                "terminal": {"state": "fully_verified", "owner": "NONE"},
                "evidence": {}, "envelope": envelope}

    good = {"service": "verify", "call_succeeded": True, "verdict": None}
    assert not list(validator.iter_errors(record(good)))
    assert list(validator.iter_errors(record({**good, "ok": True})))
    # And a verdict that is not three-state.
    for smuggled in ("true", 1, "verified"):
        assert list(validator.iter_errors(
            record({**good, "verdict": smuggled})))


def test_a_call_that_succeeded_is_not_an_entry_this_repository_may_use():
    """The decision this repository actually acts on is three-state and rests
    on the six gates, never on whether somebody's call returned."""
    succeeded = {"ok": True, "outcome": "completed"}
    with pytest.raises(ContractError):
        verdict_of(succeeded)
    entry = {"terminal": {"state": "tangent_not_verified", "owner": "INTERNAL"},
             "evidence": {n: True for n in GATES},
             "convention": {"voigt_order": ["11", "22", "33", "12", "13", "23"]}}
    answer, why = usable_for_assembly(entry)
    assert answer.is_false()
    assert "tangent_not_verified" in why


def test_every_committed_fixture_states_a_verdict_rather_than_a_call():
    """A fixture is the product of a call that succeeded. What makes it usable
    is its gates, not the fact that the call returned."""
    from contract_paths import FIXTURES
    from contract_reader import gates_of

    for path in sorted(FIXTURES.glob("*.json")):
        payload = json.loads(path.read_text())
        history = payload.get("finite_history") or {}
        assert "ok" not in history
        assert "ok" not in payload
        gates = gates_of({"evidence": history.get("evidence") or {}})
        # Every gate is three-state; none of them is a call-success.
        for name in GATES:
            assert gates[name].is_true() or gates[name].is_false() or \
                gates[name].is_not_established()
