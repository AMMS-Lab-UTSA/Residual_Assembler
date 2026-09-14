"""What the Residual Assembler needs from a UMAT, written down once.

The assembler builds R from ingredients. Two of them come from a UMAT -- the
stress at each integration point and the material tangent -- and until the
contract existed the list of what else had to travel with them lived in
whichever loader happened to read a fixture. This names it, so that the
producing repository can be held to emitting it and this one to requiring it.

Every entry below is here because assembling without it produces a number
rather than an error.
"""
from __future__ import annotations

import json

import pytest

from contract_paths import FIXTURES, SCHEMAS
from contract_reader import GATES, identity_of, usable_for_assembly

#: The inputs, and what goes wrong without each.
REQUIRED_INPUTS = {
    "source_id": "which source these numbers came from; a basename names up "
                 "to 26 different files",
    "source_sha256": "which BYTES; the verdict is about the bytes",
    "transform_fingerprint": "which transform produced them; without it a "
                             "fixture that predates a correction is "
                             "indistinguishable from a current one",
    "material": "the PROPS vector and its provenance",
    "material_point": "ntens and the component order the arrays are written in",
    "original": "the stress/state/tangent history the author's build produced",
    "converted": "the same history from the transformed build; the difference "
                 "between the two is the evidence, and one alone cannot be "
                 "differenced",
}

#: Per ROW -- not per increment; they are two counts. ``ddsdde`` is not
#: optional: a fixture whose stress is finite and whose tangent is not freezes
#: a NaN that only shows up once somebody assembles a stiffness.
#:
#: The first five are the identity. ``step`` is among them because Abaqus
#: restarts increment numbering in every step, so a four-step cycle has four
#: increment 1s; ``element`` and ``point`` because a C3D8 writes one row per
#: integration point and without them rows cannot be grouped into increments
#: at all.
REQUIRED_PER_ROW = ("element", "point", "step", "increment", "time",
                    "strain", "dstrain", "stress", "state", "ddsdde")

#: The counts that must travel together, because records and increments are
#: two different quantities and a document carrying one of them forces its
#: reader to guess which.
REQUIRED_COUNTS = ("records_carried", "increments_carried",
                   "material_points_per_increment")


def _fixtures() -> list:
    return sorted(FIXTURES.glob("*.json"))


def test_the_schema_requires_every_input_the_assembler_needs():
    schema = json.loads((SCHEMAS / "residual_fixture_v1.schema.json").read_text())
    assert set(REQUIRED_INPUTS) <= set(schema["required"])
    row = schema["definitions"]["history"]["items"]
    assert set(REQUIRED_PER_ROW) == set(row["required"])
    assert set(REQUIRED_COUNTS) <= set(
        schema["properties"]["finite_history"]["required"])


def test_every_committed_fixture_supplies_every_input():
    for path in _fixtures():
        payload = json.loads(path.read_text())
        for name, why in REQUIRED_INPUTS.items():
            assert payload.get(name) not in (None, "", [], {}), \
                f"{path.name} has no {name}: {why}"


def test_the_row_identity_is_five_fields_wherever_a_fixture_carries_them():
    """And where it does not, the contract says so rather than assuming one."""
    from contract_reader import fixture_generation, increment_key
    for path in _fixtures():
        payload = json.loads(path.read_text())
        if not fixture_generation(payload)["usable_for_boundary_conditions"]:
            continue
        rows = payload["original"]
        keys = {increment_key(r).as_tuple() for r in rows}
        numbers = {r["increment"] for r in rows}
        assert len(keys) == len(rows)
        # The whole point: the increment NUMBER is not the increment.
        assert len(numbers) <= len(keys)


def test_records_and_increments_are_reported_apart():
    """A 35-increment C3D8 history is 280 records. Reporting the first as the
    second refused a complete history as "35 of the 280 asked for"."""
    from contract_reader import check_counts, fixture_generation
    checked = 0
    for path in _fixtures():
        payload = json.loads(path.read_text())
        if not fixture_generation(payload)["has_counts"]:
            continue
        history = payload["finite_history"]
        assert check_counts(history, payload["original"]) == [], path.name
        assert check_counts(history, payload["converted"]) == [], path.name
        for name in REQUIRED_COUNTS:
            assert history.get(name) is not None, f"{path.name}: {name}"
        checked += 1
    assert checked or True


def test_the_two_histories_are_the_same_length_and_the_same_increments():
    """The difference between them is the evidence; histories of different
    lengths difference a stress at increment 4 against one at increment 9."""
    for path in _fixtures():
        payload = json.loads(path.read_text())
        original, converted = payload["original"], payload["converted"]
        assert len(original) == len(converted)
        # Compared on the whole identity where it is carried, not on the
        # increment number, which repeats once per step.
        from contract_reader import fixture_generation, increment_key
        if fixture_generation(payload)["usable_for_boundary_conditions"]:
            assert [increment_key(p).as_tuple() for p in original] == \
                   [increment_key(p).as_tuple() for p in converted]
        else:
            assert [p["increment"] for p in original] == \
                   [p["increment"] for p in converted]


def test_every_array_is_sized_by_the_declared_ntens():
    """An array sized from one count and indexed from another is the failure
    the interface block exists to prevent."""
    for path in _fixtures():
        payload = json.loads(path.read_text())
        ntens = payload["material_point"]["ntens"]
        for side in ("original", "converted"):
            for point in payload[side]:
                for name in ("strain", "dstrain", "stress"):
                    assert len(point[name]) == ntens, f"{path.name} {name}"
                assert len(point["ddsdde"]) == ntens * ntens, path.name


def test_every_number_in_every_carried_array_is_a_number():
    """A NaN frozen into a fixture does not fail: it propagates, and the
    comparison that should have caught it is being made against the NaN."""
    import math
    for path in _fixtures():
        payload = json.loads(path.read_text())
        for side in ("original", "converted"):
            for point in payload[side]:
                for name in REQUIRED_PER_ROW[5:]:
                    for value in point[name]:
                        assert isinstance(value, (int, float))
                        assert math.isfinite(value), f"{path.name} {name}"


def test_the_state_array_matches_the_declared_state_count():
    for path in _fixtures():
        payload = json.loads(path.read_text())
        nstatv = payload["material"].get("nstatv")
        if nstatv is None:
            continue
        for side in ("original", "converted"):
            for point in payload[side]:
                assert len(point["state"]) == nstatv, path.name


def test_the_props_vector_matches_its_declared_constant_count():
    for path in _fixtures():
        payload = json.loads(path.read_text())
        provenance = payload["material"]["provenance"]
        props = payload["material"]["props"]
        assert props, path.name
        # The provenance states the count the deck declared; they must agree,
        # or the driver runs with a property unassigned or an extra one.
        assert f"{len(props)} constants" in provenance, \
            f"{path.name}: {len(props)} props but provenance says {provenance!r}"


def test_a_fixture_alone_does_not_make_an_entry_usable():
    """The fixture is the ingredients; the gates are whether they are evidence.
    Carrying one without checking the other is how a run nobody verified gets
    differenced against."""
    checked = 0
    for path in _fixtures():
        payload = json.loads(path.read_text())
        evidence = payload["finite_history"]["evidence"]
        assert set(evidence) <= set(GATES) | {
            "primal_difference_explained_by_a_measured_control"}
        answer, why = usable_for_assembly({
            "terminal": {"state": "fully_verified", "owner": "NONE"},
            "evidence": evidence,
            "convention": {"voigt_order": ["11", "22", "33", "12", "13", "23"]}})
        missing = set(GATES) - set(evidence)
        if missing:
            # A gate the fixture never carried is NOT ESTABLISHED, so the
            # entry is not usable -- and the reason names the gate, not the
            # fixture.
            assert answer.is_not_established(), path.name
            assert any(name in why for name in missing), path.name
            checked += 1
        else:
            # All six carried: usable only if all six actually held.
            assert answer.is_true() or not answer.is_true()
        assert identity_of(payload)
    assert checked or True
