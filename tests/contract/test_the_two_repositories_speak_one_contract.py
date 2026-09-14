"""The shared files are the SAME files, and both ends say which version they speak.

This repository and UMAT_source_transformation exchange records. Before the
contract, each read the other's ad-hoc JSON keys: a key renamed on one side
surfaced on the other as a ``.get()`` returning ``None`` and a wrong number
carried three layers forward.

Two guards, for two different failures. The **lock** catches a shared file
edited in this repository without the version being bumped. The **handshake**
catches the two repositories being at different versions. Neither substitutes
for the other -- a repository can only digest its own copies -- and this file
holds both, plus an opportunistic byte comparison when the other checkout is
on the machine.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from contract_paths import SCHEMAS, SHARED, load_schema, umat_repo
from contract_reader import (CONTRACT_VERSION, SPEAKER, ContractVersionError,
                             handshake, require_compatible)

SHARED_SCHEMAS = ("umat_contract_v1.schema.json",
                  "residual_fixture_v1.schema.json",
                  "contract_error_v1.schema.json",
                  "transform_generation.json")


def _normalised_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def _local(name: str) -> Path:
    return (SCHEMAS / name) if name.endswith(".json") else (SHARED / name)


# ---------------------------------------------------------------------------
def test_this_repository_states_which_contract_it_speaks():
    assert handshake() == {"speaker": SPEAKER,
                           "contract_version": CONTRACT_VERSION}
    assert SPEAKER == "Residual_Assembler"


def test_every_shared_schema_declares_the_version_this_repository_speaks():
    for name in SHARED_SCHEMAS:
        document = json.loads((SCHEMAS / name).read_text())
        assert document["x-contract-version"] == CONTRACT_VERSION, name


def test_the_lock_matches_the_files_in_this_repository(lock):
    """A shared contract file edited here without the lock being regenerated
    is a change the other repository will never see, and the two ends stop
    meaning the same thing by the same key."""
    assert lock["contract_version"] == CONTRACT_VERSION
    stale = [name for name, digest in lock["files"].items()
             if _normalised_digest(_local(name)) != digest]
    assert not stale, (
        f"shared contract file(s) changed without the lock being "
        f"regenerated: {', '.join(sorted(stale))}. Either this is a contract "
        f"change -- bump the version in both repositories, regenerate the "
        f"lock, and tell the other side -- or it is an accident.")


def test_the_lock_covers_the_three_state_reader_as_well_as_the_schemas(lock):
    """The reader is part of the contract, not an implementation detail of one
    side. A reimplementation of it on this side is exactly where "null reads
    as a pass" comes back."""
    assert "tristate.py" in lock["files"]
    assert set(SHARED_SCHEMAS) < set(lock["files"])


def test_this_repository_does_not_reimplement_the_three_state_reader():
    import tristate
    assert Path(tristate.__file__).parent == SHARED
    source = (SHARED / "tristate.py").read_text()
    assert "SHARED VERBATIM" in source


def test_a_major_version_mismatch_is_refused_with_a_readable_message():
    ahead = f"{int(CONTRACT_VERSION.split('.')[0]) + 1}.0.0"
    with pytest.raises(ContractVersionError) as exc:
        require_compatible(ahead, speaker="UMAT_source_transformation")
    message = str(exc.value)
    assert ahead in message and CONTRACT_VERSION in message
    assert "UMAT_source_transformation" in message and SPEAKER in message
    assert "BREAKING" in message
    # The point of refusing: the failure mode of reading is a wrong number,
    # not a KeyError three frames deep.
    assert "wrong numbers" in message


def test_a_record_with_no_version_is_refused_rather_than_assumed_current():
    with pytest.raises(ContractVersionError) as exc:
        require_compatible(None, speaker="UMAT_source_transformation")
    assert "declares no contract version" in str(exc.value)
    with pytest.raises(ContractVersionError):
        require_compatible("1.0", speaker="UMAT_source_transformation")


def test_a_one_x_producer_is_refused_by_this_two_x_consumer():
    """The bump is not ceremonial. A 1.x producer emits rows named by two
    fields where this reader needs five, one count where it needs five, and a
    terminal enum missing the two states that were added after three entries
    which RAN were filed as runs that never happened. Every one of those a 1.x
    reader would MISREAD rather than fail on."""
    with pytest.raises(ContractVersionError) as exc:
        require_compatible("1.0.0", speaker="UMAT_source_transformation")
    assert "BREAKING" in str(exc.value)
    assert "wrong numbers" in str(exc.value)


def test_the_same_version_needs_no_caveat():
    assert require_compatible(CONTRACT_VERSION,
                              speaker="UMAT_source_transformation") == ""


# ---------------------------------------------------------------------------
def test_the_shared_files_are_byte_identical_to_the_producing_repository():
    """Opportunistic: it can only run where both checkouts are on the machine.
    When it cannot, the lock above is what holds, and this is reported as not
    established rather than passing quietly."""
    other = umat_repo()
    if other is None:
        pytest.skip(
            "UMAT_source_transformation is not on this machine, so the "
            "cross-repository byte comparison is NOT ESTABLISHED. The lock in "
            "schemas/contract_lock.json still holds this repository's copies "
            "to the digests they were published at.")
    theirs = other / "src" / "umat_oti" / "contract"
    differing = []
    for name in SHARED_SCHEMAS + ("tristate.py",):
        far = (theirs / "schemas" / name) if name.endswith(".json") \
            else (theirs / name)
        if not far.is_file():
            differing.append(f"{name} is missing from {theirs}")
        elif _normalised_digest(far) != _normalised_digest(_local(name)):
            differing.append(f"{name} differs between the two repositories")
    assert not differing, "; ".join(differing)


def test_both_repositories_record_the_same_lock():
    other = umat_repo()
    if other is None:
        pytest.skip("UMAT_source_transformation is not on this machine")
    far = json.loads((other / "src" / "umat_oti" / "contract" /
                      "contract_lock.json").read_text())
    near = json.loads((SCHEMAS / "contract_lock.json").read_text())
    assert far["combined"] == near["combined"]
    assert far["contract_version"] == near["contract_version"] == CONTRACT_VERSION


def test_the_breaking_change_policy_is_written_down():
    """A version number nobody defined the meaning of is a number."""
    policy = load_schema("umat_contract_v1")["x-breaking-change"].lower()
    for rule in ("removing", "renaming", "required", "meaning"):
        assert rule in policy
