"""A fixture whose transform is not the store's is evidence about a transform
that no longer exists.

A regression fixture is what later runs are compared against, so anything
wrong inside it is wrong in every comparison made against it afterwards,
silently. The way it goes wrong here is not a NaN -- the exporter refuses
those. It is a fixture that was correct when it was frozen and stopped being
evidence when the transform changed underneath it, and which looks exactly
like a current one. The frozen collection was found to hold 67 materials whose
tangent numbers predated a correction.

So every fixture carries the transform fingerprint it was cut at, and this
repository -- the consumer, the one that would difference against it -- is
able to refuse one that does not match. Not the exporter alone: a claim
checked only by the tool that makes it is a claim nobody checked.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from contract_paths import CURRENT_FIXTURES, FIXTURES, SCHEMAS
from contract_reader import (CONTRACT_VERSION, ContractError,
                             check_fingerprint, fixture_generation,
                             identity_of, require_current,
                             require_five_field_identity)
from tristate import Tri

#: The transform this contract's evidence belongs to, read from the shared
#: generation file rather than written down here as well. A constant
#: hard-coded on each side of the boundary is two constants, and the day they
#: drift this repository accepts a fixture the producer would have refused.
STORE_FINGERPRINT = json.loads(
    (SCHEMAS / "transform_generation.json").read_text())["transform_fingerprint"]

#: A fingerprint that is definitely not the current one, for exercising the
#: rejection path whichever fixture set happens to be checked out. It is the
#: real generation the first four fixtures here were frozen at.
A_SUPERSEDED_FINGERPRINT = "ff94800b1884bcc0"

FIXTURE_SCHEMA = "umat-oti/residual-fixture/1"


def _fixtures() -> list:
    return sorted(FIXTURES.glob("*.json"))


def _validate(payload: dict, schema_name: str, where: str) -> None:
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((SCHEMAS / f"{schema_name}.schema.json").read_text())
    errors = sorted(jsonschema.Draft7Validator(schema).iter_errors(payload),
                    key=lambda e: list(e.path))
    assert not errors, f"{where}: " + "; ".join(
        f"at {'/'.join(str(p) for p in e.absolute_path) or '(root)'}: "
        f"{e.message}" for e in errors[:5])


# ---------------------------------------------------------------------------
def test_a_matching_fingerprint_is_the_only_pass():
    assert check_fingerprint(STORE_FINGERPRINT, STORE_FINGERPRINT).is_true()
    stale = check_fingerprint(A_SUPERSEDED_FINGERPRINT, STORE_FINGERPRINT,
                              where="a.json")
    assert stale.is_false()
    assert "has since changed" in stale.why


def test_a_fixture_with_no_fingerprint_is_not_current_by_default():
    """Three states. It might be current and it might predate a correction,
    and nothing on disk says which."""
    answer = check_fingerprint(None, STORE_FINGERPRINT, where="a.json")
    assert answer.is_not_established()
    assert not answer.is_true() and not answer.is_false()
    assert "not current by default" in answer.why


def test_require_current_refuses_not_established_as_well_as_stale():
    """"Nothing said" is not known-current, and treating it as one is the
    null-reads-as-a-pass defect wearing a different hat."""
    require_current(STORE_FINGERPRINT, STORE_FINGERPRINT)
    for refused in (A_SUPERSEDED_FINGERPRINT, "", None):
        with pytest.raises(ContractError):
            require_current(refused, STORE_FINGERPRINT)


# ---------------------------------------------------------------------------
# what is actually committed here
# ---------------------------------------------------------------------------
def test_the_generation_file_is_the_single_source_of_the_fingerprint():
    """Both repositories read it; neither writes it down twice."""
    generation = json.loads((SCHEMAS / "transform_generation.json").read_text())
    assert generation["transform_fingerprint"] == STORE_FINGERPRINT
    assert generation["x-contract-version"] == CONTRACT_VERSION
    assert "NOT ESTABLISHED" in generation["consumers_must"]
    assert STORE_FINGERPRINT != A_SUPERSEDED_FINGERPRINT


def test_some_fixtures_are_committed():
    assert _fixtures()


def test_current_operational_fixtures_carry_independent_verification():
    from residual_core.materials.verified_fixture import load

    current = [path for path in sorted(CURRENT_FIXTURES.glob("*.json"))
               if json.loads(path.read_text())["transform_fingerprint"]
               == STORE_FINGERPRINT]
    assert current, "regenerate operational fixtures with the current producer"
    assert {path.name for path in current} == {
        "isotropic-elasticity--f7eb90376a.json", "j2_props--2feae9f158.json"}
    for path in current:
        payload = json.loads(path.read_text())
        archived = json.loads((FIXTURES / path.name).read_text())
        fixture = load(path)
        assert fixture.all_six_gates, path.name
        require_five_field_identity(payload, where=path.name)
        _validate(payload, "residual_fixture_v1", path.name)
        for field in ("source_id", "source_sha256", "deck", "deck_digest", "material"):
            assert payload[field] == archived[field], (path.name, field)
        assert len(payload["original"]) == len(archived["original"])
        assert fixture.verification["states_checked"] > 0
        assert fixture.verification["states_agreeing"] > 0
        for side in ("original", "transformed"):
            scan = payload["finite_history"]["whole_history"][side]
            assert scan["values_scanned"] > 0
            assert scan["first_non_finite"] is None
            grouping = payload["finite_history"]["history_grouping"][side]
            assert grouping == archived["finite_history"]["history_grouping"][side]


def test_historical_archive_preserves_every_original_inventory_hash():
    inventory = json.loads((SCHEMAS.parent / "docs/evidence/recovery_evidence_inventory.json").read_text())
    artifacts = inventory["repositories"]["RA"]["artifacts"]
    assert len(artifacts) == len(_fixtures()) == 10
    for artifact in artifacts:
        path = FIXTURES / Path(artifact["path"]).name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == artifact["sha256"]
        assert path.stat().st_size == artifact["bytes"]


@pytest.mark.parametrize("name,held,unknown", [
    ("isotropic-elasticity--f7eb90376a.json", 1, 2),
    ("j2_props--2feae9f158.json", 32, 3),
])
def test_current_fixture_cli_uses_default_currency_and_real_references(
        tmp_path, name, held, unknown):
    from residual_core.core.fixture_residual_check import check_fixture

    output = tmp_path / "report.json"
    result = subprocess.run(
        [sys.executable, "-m", "residual_core.core.fixture_residual_check",
         "--fixture", str(CURRENT_FIXTURES / name), "--out", str(output)],
        cwd=SCHEMAS.parent, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads(output.read_text())
    assert not report["refused"]
    assert "fingerprint_note" not in report
    assert report["fixture_fingerprint"] == STORE_FINGERPRINT
    assert report["summary"]["held"] == held
    assert report["summary"]["failed"] == 0
    assert report["summary"]["not_established"] == unknown
    if name.startswith("j2_props"):
        assert "umat_oti.validation.j2_reference" in report["material_model"]
        assert {check["name"] for check in report["checks"]
                if check["status"] == "holds"} >= {
                    "dR/du", "dR/dq", "dR/dE", "dR/dNU", "dR/dSIGY0", "dR/dH"}
    refused = check_fixture(FIXTURES / name)
    assert "frozen under transform fingerprint" in refused["refused"]
    assert refused["checks"] == []


def test_retired_pass12_fixtures_are_history_not_current_baselines():
    from residual_core.materials.verified_fixture import (
        CURRENT_TRANSFORM_FINGERPRINT, STORE_PROVENANCE, FixtureError, load,
    )

    assert CURRENT_TRANSFORM_FINGERPRINT == STORE_FINGERPRINT
    assert STORE_PROVENANCE["fingerprint"] == "94a92c01814f107a"
    assert STORE_PROVENANCE["status"] == "historical_only"
    assert STORE_PROVENANCE["fingerprint"] != CURRENT_TRANSFORM_FINGERPRINT
    retired = []
    for path in _fixtures():
        payload = json.loads(path.read_text())
        if payload["transform_fingerprint"] != "94a92c01814f107a":
            continue
        retired.append(path.name)
        assert payload["original"] and payload["converted"]
        identity_of(payload)
        with pytest.raises(ContractError) as exc:
            require_current(payload["transform_fingerprint"], STORE_FINGERPRINT,
                            where=path.name)
        assert path.name in str(exc.value)
        assert payload["transform_fingerprint"] in str(exc.value)
        assert STORE_FINGERPRINT in str(exc.value)
        with pytest.raises(FixtureError) as refused:
            load(path)
        assert path.name in str(refused.value)
        assert payload["transform_fingerprint"] in str(refused.value)
        assert STORE_FINGERPRINT in str(refused.value)
    assert retired, "the pass12 evidence must remain available as history"


def test_a_two_x_fixture_validates_and_a_one_x_one_says_why_it_does_not():
    """The contract failing clearly, rather than a consumer discovering the
    missing step three layers away.

    Measured: one committed fixture speaks 2.x -- five identity fields and
    five counts -- and the rest speak 1.x. A 1.x fixture is not broken and is
    not current either.
    """
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((SCHEMAS / "residual_fixture_v1.schema.json").read_text())
    validator = jsonschema.Draft7Validator(schema)
    modern, legacy = [], []
    for path in _fixtures():
        payload = json.loads(path.read_text())
        generation = fixture_generation(payload)
        errors = list(validator.iter_errors(payload))
        if generation["usable_for_boundary_conditions"]:
            modern.append(path.name)
            assert not errors, f"{path.name}: {errors[0].message}"
            require_five_field_identity(payload, where=path.name)
        else:
            legacy.append(path.name)
            assert errors, f"{path.name} lacks the five fields and validated"
            assert any(f in e.message for e in errors
                       for f in ("element", "point", "step"))
            with pytest.raises(ContractError) as exc:
                require_five_field_identity(payload, where=path.name)
            assert "3.0 relative" in str(exc.value)
            # And its numbers are still readable -- refusing a usable fixture
            # is its own kind of wrong answer.
            assert payload["original"] and payload["converted"]
    assert modern or legacy
    assert len(modern) + len(legacy) == len(_fixtures())


def test_every_committed_fixture_is_identified_by_path_and_digest():
    """Not by basename. 23 basenames in the frozen corpus name more than one
    source file, and two of these four fixtures are called ``umat``."""
    identities = [identity_of(json.loads(p.read_text())) for p in _fixtures()]
    assert identities
    for path, digest in identities:
        assert "/" in path
        assert len(digest) == 64
    # And they are distinct, which a basename would not have been: two of the
    # committed fixtures are called "umat".
    assert len({p for p, _ in identities}) == len(identities)
    assert len({d for _, d in identities}) == len(identities)


def test_every_committed_fixture_declares_the_transform_it_was_frozen_at():
    for path in _fixtures():
        payload = json.loads(path.read_text())
        assert payload.get("transform_fingerprint"), path.name
        assert payload["schema"] == FIXTURE_SCHEMA


def test_every_committed_fixture_is_either_current_or_refused_by_name():
    """The rule, held against whatever is actually committed here.

    Deliberately not a hard-coded count or a hard-coded fingerprint. Two
    fixture sets exist in this project and they are in different states: four
    frozen at ff94800b1884bcc0, which this contract refuses as regression
    baselines, and nine since re-frozen at the current generation, which it
    accepts. A test asserting either state would start failing the moment the
    other was checked out, and what is worth holding is the rule -- every
    fixture says which transform produced it, a stale one is refused with both
    fingerprints named, and a current one is accepted.
    """
    current, stale, unfingerprinted = [], [], []
    for path in _fixtures():
        payload = json.loads(path.read_text())
        their_fingerprint = payload.get("transform_fingerprint")
        answer = check_fingerprint(their_fingerprint, STORE_FINGERPRINT,
                                   where=path.name)
        if answer.is_true():
            current.append(path.name)
            require_current(their_fingerprint, STORE_FINGERPRINT,
                            where=path.name)
        elif answer.is_false():
            stale.append(path.name)
            with pytest.raises(ContractError) as exc:
                require_current(their_fingerprint, STORE_FINGERPRINT,
                                where=path.name)
            # Both fingerprints named, so a reader knows what to re-freeze at.
            assert str(their_fingerprint) in str(exc.value)
            assert STORE_FINGERPRINT in str(exc.value)
        else:
            unfingerprinted.append(path.name)

    # Exhaustive: a fixture carrying no fingerprint cannot be told from one
    # that predates a correction, and shipping one is shipping a baseline
    # nobody can date.
    assert not unfingerprinted, (
        f"{len(unfingerprinted)} committed fixture(s) carry no "
        f"transform_fingerprint: {', '.join(unfingerprinted)}")
    assert len(current) + len(stale) == len(_fixtures())


def test_a_superseded_fixture_is_refused_even_when_everything_else_is_fine():
    """The rejection path, exercised regardless of which set is checked out."""
    payload = json.loads(_fixtures()[0].read_text())
    identity_of(payload)                       # everything else about it is fine
    answer = check_fingerprint(A_SUPERSEDED_FINGERPRINT, STORE_FINGERPRINT,
                               where="a.json")
    assert answer.is_false()
    assert "evidence about the old transform" in answer.why
    with pytest.raises(ContractError):
        require_current(A_SUPERSEDED_FINGERPRINT, STORE_FINGERPRINT)


def test_a_superseded_fixture_is_still_readable_as_history():
    """Refusing one as a REGRESSION BASELINE is not the same as deleting it.
    These fixtures are the only end-to-end evidence that the assembler
    consumes what the pipeline emits, and the contract's answer to a stale one
    is "not current", not "not a fixture"."""
    for path in _fixtures():
        payload = json.loads(path.read_text())
        assert payload["original"] and payload["converted"]
        assert len(payload["original"]) == len(payload["converted"])
        identity_of(payload)          # still a legal identity


def test_a_fixture_carrying_a_different_fingerprint_from_its_own_store_is_refused(tmp_path):
    """The general rule, on a fixture built for the purpose rather than on the
    four that happen to be stale."""
    payload = json.loads(_fixtures()[0].read_text())
    payload["transform_fingerprint"] = STORE_FINGERPRINT
    assert check_fingerprint(payload["transform_fingerprint"],
                             STORE_FINGERPRINT).is_true()
    payload["transform_fingerprint"] = "0123456789abcdef"
    with pytest.raises(ContractError) as exc:
        require_current(payload["transform_fingerprint"], STORE_FINGERPRINT,
                        where="rebuilt.json")
    assert "0123456789abcdef" in str(exc.value)


def test_a_fixture_with_no_fingerprint_field_at_all_fails_the_schema(tmp_path):
    jsonschema = pytest.importorskip("jsonschema")
    payload = json.loads(_fixtures()[0].read_text())
    del payload["transform_fingerprint"]
    schema = json.loads((SCHEMAS / "residual_fixture_v1.schema.json").read_text())
    errors = list(jsonschema.Draft7Validator(schema).iter_errors(payload))
    assert errors
    assert any("transform_fingerprint" in e.message for e in errors)


def test_a_fixture_identified_by_a_basename_fails_the_schema():
    jsonschema = pytest.importorskip("jsonschema")
    payload = json.loads(_fixtures()[0].read_text())
    payload["source_id"] = "umat.f"
    schema = json.loads((SCHEMAS / "residual_fixture_v1.schema.json").read_text())
    assert list(jsonschema.Draft7Validator(schema).iter_errors(payload))
    with pytest.raises(ContractError) as exc:
        identity_of(payload)
    assert "basename" in str(exc.value)


def test_a_material_in_a_fixture_never_travels_without_provenance():
    """A fixture is compared against forever after; a constant in one with no
    provenance is indistinguishable from an invented one."""
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((SCHEMAS / "residual_fixture_v1.schema.json").read_text())
    for path in _fixtures():
        payload = json.loads(path.read_text())
        provenance = payload["material"]["provenance"]
        assert provenance.strip(), path.name
        assert "*MATERIAL" in provenance
        payload["material"]["provenance"] = ""
        assert list(jsonschema.Draft7Validator(schema).iter_errors(payload))


def test_the_tensor_convention_is_stated_in_every_fixture():
    """Engineering shear read as tensorial is wrong by a factor of two in
    three components and right everywhere else."""
    for path in _fixtures():
        point = json.loads(path.read_text())["material_point"]
        assert point["ntens"] >= 1
        assert "ENGINEERING" in point["voigt_order"].upper()
