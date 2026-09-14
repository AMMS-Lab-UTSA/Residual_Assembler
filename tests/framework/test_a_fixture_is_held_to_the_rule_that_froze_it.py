"""A fixture is what later runs are compared against, so the loader has to
hold it to the same rule the exporter does -- and to the SAME rule, not to an
older copy of it.

``export_residual_fixture.py`` in UMAT_source_transformation refuses to write a
fixture that is not a completely finite successful history. ``load`` here
refuses to read one, because a claim checked only by the tool that makes it is
a claim nobody checked. Two copies of one rule are two rules and they drift,
and they had: the exporter grew two refusals this loader did not have.

**A short window.** The exporter's window clamps itself to whatever is on disk,
so an analysis that stopped at three increments produced a three-increment
fixture. Nothing comparing against it could tell that from a run meant to be
that length. The exporter now writes ``increments_requested`` beside
``increments_carried``; this reads both.

**A NaN outside the carried window.** Six increments out of hundreds are finite
in a run that went to NaN at increment 200, so scanning only what is carried
cannot tell a verification from the salvageable prefix of a failed analysis.
The exporter now scans the whole history and writes what it found; this reads
that rather than taking the window's finiteness as evidence about the run.

And the third thing, which is the reason this file is not simply stricter: the
four fixtures committed here were frozen before either check existed. They do
not carry the keys. Refusing them would be treating "frozen before the check"
as "failed the check", and passing them silently would be treating it as
"passed". So the loader names which claims a fixture carries and which it does
not, and the two are different answers.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from residual_core.materials.verified_fixture import (CHECKABLE_CLAIMS,
                                                      FixtureError, load,
                                                      load_all)

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "verified"


@pytest.fixture()
def payload() -> dict:
    return json.loads((FIXTURES / "j2_props--2feae9f158.json").read_text())


def _written(tmp_path: Path, payload: dict, name: str = "case.json") -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# what is committed here
# ---------------------------------------------------------------------------
def test_every_committed_fixture_still_loads():
    """Measured: four fixtures, each six increments of original and converted,
    every number a number."""
    fixtures = load_all(FIXTURES)
    assert len(fixtures) == 4
    for fixture in fixtures:
        assert fixture.increments() == 6
        assert fixture.ntens > 0
        assert fixture.source_id


def test_a_fixture_says_which_claims_it_carries_and_which_it_does_not():
    """"Checked and clean" and "never checked" are different answers, and a
    consumer that cannot tell them apart is treating an absence as a pass.

    The four committed here predate the two newest checks, so each names them
    as not carried rather than reporting them as passed.
    """
    for fixture in load_all(FIXTURES):
        assert "the run was finite from end to end" in fixture.claims_checked
        assert fixture.claims_not_carried, fixture.source_id
        assert ("the window is the length the export asked for, not what was "
                "on disk") in fixture.claims_not_carried
        assert ("every number the run wrote is finite, not only the ones "
                "carried") in fixture.claims_not_carried
        assert len(fixture.claims_checked) + len(fixture.claims_not_carried) \
            == len(CHECKABLE_CLAIMS) == 5


def test_a_fixture_carrying_every_claim_reports_none_missing(tmp_path, payload):
    payload["finite_history"].update({
        "increments_requested": 6, "records_available": 280,
        "whole_history": {"original": {"first_non_finite": None,
                                       "records_scanned": 280},
                          "transformed": {"first_non_finite": None,
                                          "records_scanned": 280}}})
    fixture = load(_written(tmp_path, payload))
    assert fixture.claims_not_carried == ()
    assert len(fixture.claims_checked) == len(CHECKABLE_CLAIMS)


# ---------------------------------------------------------------------------
# the short window
# ---------------------------------------------------------------------------
def test_a_window_shorter_than_the_export_asked_for_is_refused(tmp_path,
                                                               payload):
    """The refusal the exporter grew. A three-increment fixture out of a
    six-increment request is a fixture built on however far an analysis got."""
    payload["finite_history"].update({"increments_requested": 6,
                                      "increments_carried": 3,
                                      "records_available": 3})
    payload["original"] = payload["original"][:3]
    payload["converted"] = payload["converted"][:3]
    with pytest.raises(FixtureError) as raised:
        load(_written(tmp_path, payload))
    said = str(raised.value)
    assert "3 increment(s) of the 6" in said
    assert "however far an analysis got before it stopped" in said


def test_a_window_of_exactly_what_was_asked_for_is_not_refused(tmp_path,
                                                               payload):
    """The complaint is about a silent shortfall, not about a small fixture."""
    payload["finite_history"].update({"increments_requested": 6,
                                      "increments_carried": 6,
                                      "records_available": 280})
    assert load(_written(tmp_path, payload)).increments() == 6


def test_a_fixture_that_does_not_state_what_was_asked_for_is_not_refused(
        tmp_path, payload):
    """It is reported as not carrying the claim, which is the third answer."""
    payload["finite_history"].pop("increments_requested", None)
    fixture = load(_written(tmp_path, payload))
    assert "the window is the length the export asked for, not what was on " \
           "disk" in fixture.claims_not_carried


# ---------------------------------------------------------------------------
# the NaN the window cannot see
# ---------------------------------------------------------------------------
def test_a_run_that_went_non_finite_outside_the_window_is_refused(tmp_path,
                                                                  payload):
    """Every number in this fixture is finite. The run it came out of was not,
    and only the exporter's whole-history scan knows that."""
    payload["finite_history"]["whole_history"] = {
        "original": {"records_scanned": 280, "values_scanned": 11480,
                     "first_non_finite": {"record": 201, "array": "STRESS",
                                          "index": 3, "value": "nan",
                                          "increment": 26, "element": 1,
                                          "point": 5}},
        "transformed": {"records_scanned": 280, "first_non_finite": None}}
    with pytest.raises(FixtureError) as raised:
        load(_written(tmp_path, payload))
    said = str(raised.value)
    assert "OUTSIDE the carried window" in said
    assert "STRESS[3] is nan at record 201" in said
    assert "prefix of a failed analysis" in said


def test_a_non_finite_tangent_outside_the_window_is_refused_the_same_way(
        tmp_path, payload):
    payload["finite_history"]["whole_history"] = {
        "transformed": {"records_scanned": 280,
                        "first_non_finite": {"record": 90, "array": "DDSDDE",
                                             "index": 7, "value": "inf",
                                             "increment": 12}}}
    with pytest.raises(FixtureError) as raised:
        load(_written(tmp_path, payload))
    assert "DDSDDE[7] is inf" in str(raised.value)
    assert "transformed run" in str(raised.value)


def test_a_clean_whole_history_scan_is_not_refused(tmp_path, payload):
    payload["finite_history"]["whole_history"] = {
        "original": {"records_scanned": 280, "first_non_finite": None},
        "transformed": {"records_scanned": 280, "first_non_finite": None}}
    assert load(_written(tmp_path, payload)).increments() == 6


# ---------------------------------------------------------------------------
# the refusals that were already here still fire
# ---------------------------------------------------------------------------
def test_a_nan_among_the_carried_numbers_is_still_refused(tmp_path, payload):
    payload["original"][2]["stress"][1] = float("nan")
    with pytest.raises(FixtureError) as raised:
        load(_written(tmp_path, payload))
    assert "stress[1]" in str(raised.value)


def test_a_run_that_was_not_finite_end_to_end_is_still_refused(tmp_path,
                                                               payload):
    payload["finite_history"]["complete_finite_verification_run"] = False
    with pytest.raises(FixtureError) as raised:
        load(_written(tmp_path, payload))
    assert "complete_finite_verification_run" in str(raised.value)


def test_a_gate_that_was_never_measured_is_not_taken_as_passed(tmp_path,
                                                               payload):
    payload["finite_history"]["evidence"].pop("complete_history_finite")
    with pytest.raises(FixtureError) as raised:
        load(_written(tmp_path, payload))
    assert "complete_history_finite='not measured'" in str(raised.value)
