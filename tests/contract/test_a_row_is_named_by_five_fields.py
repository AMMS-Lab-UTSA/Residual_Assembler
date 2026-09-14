"""Five fields name a row here too, and records are not increments.

Both halves are checked on the consuming side and not only where the fixture
was written, because a claim checked only by the tool that makes it is a claim
nobody checked -- and this repository is the one that would be wrong. It
rebuilds boundary conditions from a fixture's rows and it integrates over the
points of an element, so it is the side that pays for either mistake.

The identity module is imported from ``_shared``, a verbatim copy of the
producer's, rather than reimplemented: two implementations of one rule are two
rules, and the last time these two ends disagreed about which count was which,
a complete 280-record history was refused as "35 of the 280 asked for".
"""
from __future__ import annotations

import json

import pytest

from contract_paths import FIXTURES
from contract_reader import (COUNT_FIELDS, IDENTITY_FIELDS, ContractError,
                             FrameError, check_counts, count_history,
                             fixture_generation, frame_key, group_by_increment,
                             increment_key, point_key,
                             require_five_field_identity)


def _fixtures() -> list:
    return sorted(FIXTURES.glob("*.json"))


def _five_field():
    for path in _fixtures():
        payload = json.loads(path.read_text())
        if fixture_generation(payload)["usable_for_boundary_conditions"]:
            return path, payload
    return None, None


# ---------------------------------------------------------------------------
def test_the_identity_is_five_fields():
    assert IDENTITY_FIELDS == ("element", "point", "step", "increment", "time")


def test_a_row_missing_any_one_is_refused_by_name():
    complete = {"element": 1, "point": 1, "step": 1, "increment": 1, "time": 0.0}
    frame_key(complete)
    for dropped in IDENTITY_FIELDS:
        with pytest.raises(FrameError) as exc:
            frame_key({k: v for k, v in complete.items() if k != dropped})
        assert dropped in str(exc.value)


def test_four_increment_ones_are_four_different_rows():
    """Abaqus restarts increment numbering in every step. A four-step cycle
    therefore has four increment 1s, and they are four different states of the
    material -- reconstructing the loading from the wrong one was wrong by 3.0
    relative, a different deformation and not a tolerance."""
    rows = [{"element": 1, "point": 1, "step": s, "increment": 1,
             "time": float(s - 1)} for s in (1, 2, 3, 4)]
    assert len({r["increment"] for r in rows}) == 1
    assert len({increment_key(r).as_tuple() for r in rows}) == 4
    assert len(group_by_increment(rows)) == 4


def test_a_point_is_not_an_increment():
    """Mixing them makes two integration points of one increment look like two
    increments, which is how a 280-record history became 280 increments."""
    rows = [{"element": 1, "point": p, "step": 1, "increment": 1, "time": 0.1}
            for p in range(1, 9)]
    assert len(group_by_increment(rows)) == 1
    assert len({point_key(r).as_tuple() for r in rows}) == 8
    counts = count_history(rows)
    assert counts.records_carried == 8
    assert counts.increments_carried == 1


def test_a_c3d8_history_is_more_records_than_increments():
    rows = [{"element": 1, "point": p, "step": 1, "increment": i,
             "time": i * 0.1} for i in range(1, 36) for p in range(1, 9)]
    counts = count_history(rows, material_points_available=8)
    assert (counts.records_carried, counts.increments_carried) == (280, 35)
    assert counts.material_points_per_increment == 8


def test_counts_that_contradict_each_other_are_refused():
    assert check_counts({"records_carried": 280, "increments_carried": 280,
                         "material_points_per_increment": 8})
    assert check_counts({"records_carried": 6, "increments_carried": 35,
                         "material_points_per_increment": 0})
    assert check_counts({"increments_carried": 35})      # only one count
    assert check_counts({"records_carried": 280, "increments_carried": 35,
                         "material_points_per_increment": 8}) == []


def test_the_five_counts_are_named():
    assert COUNT_FIELDS == ("records_carried", "increments_carried",
                            "material_points_per_increment",
                            "material_point_carried",
                            "material_points_available")


# ---------------------------------------------------------------------------
# what is actually committed here
# ---------------------------------------------------------------------------
def test_a_committed_five_field_fixture_shows_the_collision():
    path, payload = _five_field()
    if payload is None:
        pytest.skip("no five-field fixture committed here yet")
    rows = payload["original"]
    keys = {increment_key(r).as_tuple() for r in rows}
    numbers = {r["increment"] for r in rows}
    assert len(keys) == len(rows)
    assert len(numbers) < len(keys), (
        f"{path.name}: {len(numbers)} increment numbers for {len(rows)} rows")
    ones = [r for r in rows if r["increment"] == 1]
    assert len({r["step"] for r in ones}) == len(ones)
    assert len({round(r["stress"][0], 6) for r in ones}) == len(ones)


def test_a_committed_fixtures_counts_agree_with_its_own_rows():
    for path in _fixtures():
        payload = json.loads(path.read_text())
        if not fixture_generation(payload)["has_counts"]:
            continue
        assert check_counts(payload["finite_history"],
                            payload["original"]) == [], path.name


def test_a_one_x_fixture_is_refused_for_the_loading_and_kept_for_its_numbers():
    for path in _fixtures():
        payload = json.loads(path.read_text())
        generation = fixture_generation(payload)
        if generation["usable_for_boundary_conditions"]:
            require_five_field_identity(payload, where=path.name)
        else:
            with pytest.raises(ContractError) as exc:
                require_five_field_identity(payload, where=path.name)
            assert "3.0 relative" in str(exc.value)
            assert "may still be differenced against" in str(exc.value)
        assert payload["original"] and payload["converted"]


def test_this_repository_does_not_reimplement_the_identity_rule():
    import frames
    from pathlib import Path as _P
    assert _P(frames.__file__).parent.name == "_shared"
    assert "SHARED" in _P(frames.__file__).read_text() or \
        "shared" in frames.__doc__.lower()
