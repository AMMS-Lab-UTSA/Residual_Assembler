"""Every terminal state that reaches this repository says whose move it is.

EXTERNAL is a fact about somebody's published repository: nobody wrote down
what the material is made of, the file is not a UMAT, the source does not
compile. More engineering in the producing project changes none of them.

INTERNAL is a limitation of that project: the transform refused, the deck
generator had no element, the finite difference could not resolve one. Every
one of these is work.

Why this repository cares. The assembler's README and status reporting say how
much of the corpus is available to drive it. A count that pooled the two would
be a claim about the corpus made out of facts about the pipeline -- it would
report somebody else's file as at fault for a gap that is the producing
project's, or hide that project's work behind a published-source excuse. So
the owner travels with the state, and a state with no owner is refused rather
than defaulted.
"""
from __future__ import annotations

import json

import pytest

from contract_paths import SCHEMAS
from contract_reader import (EXTERNAL_STATES, FULLY_VERIFIED, INTERNAL_STATES,
                             ContractError, owner_of, usable_for_assembly)

#: The two states added after the silent default did its damage. Three of the
#: 237 frozen entries carry the first: they RAN, and the vocabulary's old
#: fall-through answered ``not_attempted`` -- INTERNAL, "the run never
#: happened" -- about all three.
ADDED_AFTER_THE_DEFAULT_DID_DAMAGE = {
    "arguments_diverged_before_the_routine": "INTERNAL",
    "disagreement_not_in_any_recorded_call": "INTERNAL",
}
#: A rung nothing has ever mapped, for exercising the refusal itself.
UNMAPPED_IN_THE_FROZEN_STORE = "a_rung_nobody_has_ever_mapped"


def test_the_three_owners_are_distinct_and_nothing_is_in_two():
    assert not (EXTERNAL_STATES & INTERNAL_STATES)
    assert FULLY_VERIFIED not in EXTERNAL_STATES
    assert FULLY_VERIFIED not in INTERNAL_STATES
    assert owner_of(FULLY_VERIFIED) == "NONE"
    for state in EXTERNAL_STATES:
        assert owner_of(state) == "EXTERNAL"
    for state in INTERNAL_STATES:
        assert owner_of(state) == "INTERNAL"


def test_the_two_states_added_after_the_default_did_damage_are_known_here():
    """A consumer that falls through on them is the consumer the old default
    was written for."""
    for state, owner in ADDED_AFTER_THE_DEFAULT_DID_DAMAGE.items():
        assert owner_of(state) == owner
    assert "arguments_diverged_before_the_routine" in INTERNAL_STATES
    assert "disagreement_not_in_any_recorded_call" in INTERNAL_STATES


def test_a_state_with_no_owner_is_refused_rather_than_defaulted():
    with pytest.raises(ContractError) as exc:
        owner_of(UNMAPPED_IN_THE_FROZEN_STORE)
    message = str(exc.value)
    assert "EXTERNAL" in message and "INTERNAL" in message
    assert "guessing would assert an owner nobody established" in message
    with pytest.raises(ContractError):
        owner_of("something_nobody_mapped")
    with pytest.raises(ContractError):
        owner_of("")


def test_a_pause_waiting_on_terminal_input_belongs_to_the_author():
    """It is a property of the file somebody published, not of this pipeline,
    even though what it does is hold an Abaqus licence until a timeout."""
    assert owner_of("waits_for_input") == "EXTERNAL"


def test_a_transform_refusal_is_this_projects_work_and_not_the_authors():
    assert owner_of("transform_refused") == "INTERNAL"
    assert owner_of("tangent_not_verified") == "INTERNAL"
    assert owner_of("experiment_not_generated") == "INTERNAL"


def test_a_published_stub_is_external_and_its_own_state():
    """A file presenting the 37-argument UMAT header that assigns neither
    STRESS nor DDSDDE is not a model this project failed to convert -- there
    is no model. It is not ``transform_refused`` (internal) and it is not
    ``incomplete_or_corrupt_source``, which is glossed 'does not compile',
    and a template compiles perfectly well."""
    assert owner_of("published_stub_no_constitutive_content") == "EXTERNAL"
    assert owner_of("transform_refused") == "INTERNAL"


def test_the_schema_refuses_a_state_declared_with_the_wrong_owner():
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((SCHEMAS / "umat_contract_v1.schema.json").read_text())
    validator = jsonschema.Draft7Validator(schema)

    def record(state, owner):
        return {"contract_version": "3.0.0",
                "identity": {"path": "repo__x/u.f", "sha256": "a" * 64},
                "transform_fingerprint": "b0d27ee53c630500",
                "terminal": {"state": state, "owner": owner},
                "evidence": {}}

    assert not list(validator.iter_errors(record("not_a_umat", "EXTERNAL")))
    assert not list(validator.iter_errors(record("transform_refused", "INTERNAL")))
    # An internal limitation transmitted as an external blocker, and the reverse.
    assert list(validator.iter_errors(record("transform_refused", "EXTERNAL")))
    assert list(validator.iter_errors(record("not_a_umat", "INTERNAL")))
    assert list(validator.iter_errors(record("fully_verified", "EXTERNAL")))
    # And a state the vocabulary does not have at all.
    assert list(validator.iter_errors(
        record(UNMAPPED_IN_THE_FROZEN_STORE, "INTERNAL")))
    # And the two added states validate on the side the evidence puts them.
    assert not list(validator.iter_errors(
        record("arguments_diverged_before_the_routine", "INTERNAL")))
    assert list(validator.iter_errors(
        record("arguments_diverged_before_the_routine", "EXTERNAL")))
    assert not list(validator.iter_errors(
        record("disagreement_not_in_any_recorded_call", "INTERNAL")))


def test_only_a_fully_verified_entry_may_drive_the_assembler():
    for state in sorted(EXTERNAL_STATES | INTERNAL_STATES):
        answer, why = usable_for_assembly({
            "terminal": {"state": state, "owner": owner_of(state)},
            "evidence": {}})
        assert answer.is_false(), state
        assert state in why and owner_of(state) in why
