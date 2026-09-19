"""Two implementations of one contract, asked the same questions about the
same 237 records, and required to give the same answers.

This is the test the whole contract exists to make possible. Everything else
here checks that this repository reads its own copy of the rules correctly.
This checks that reading them correctly produces the SAME ANSWER as the
producing repository gets -- which is the only property that matters, because
a contract two sides each implement correctly and differently is not a
contract.

It runs when the producing checkout is beside this one and the frozen store is
on the machine. When it cannot run it says so rather than passing quietly:
"the two readers were never compared" is a third answer and not a pass, which
is the same rule this contract applies to every gate it carries.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from contract_paths import runs_the_checkouts_code, umat_repo
from contract_reader import (GATES, SEVENTH, ContractError, gates_of,
                             identity_of, owner_of, primal_settled,
                             require_compatible, usable_for_assembly)

STORE_RELATIVE = Path("corpus_run") / "pass11" / "results" / \
    "store_verification.jsonl"


def _producer():
    """The producing checkout, or a skip that says why it was not compared.

    The import trap: whichever ``umat_oti`` is imported first is the one every
    later import gets, and other tests in this repository's suite import it
    before this module's tests run. That copy may be an editable install of
    ANOTHER checkout (the main one, not the worktree under test), or an
    install of another version altogether.

    So this refuses to compare against a producer it did not choose -- but
    "chose" means the CODE, not the directory it was imported from. A copy
    installed from this very checkout (``pip install <checkout>``, as CI and
    users install it) lives in site-packages and runs exactly the checkout's
    code; when every Python file is byte-identical to the checkout's
    ``src/umat_oti`` the comparison runs. Anything else is a different
    version, and the module skips rather than fails, because "another test
    bound another package first" is a fact about the run and not about the
    contract; and the skip says NOT COMPARED rather than passing quietly,
    which is the same rule this contract applies to every gate it carries.
    Run this module on its own, or set ``UMAT_OTI_REPO`` to the checkout the
    imported package came from, to make the comparison actually happen.
    """
    repo = umat_repo()
    if repo is None:
        pytest.skip(
            "no UMAT_source_transformation checkout beside this one, so the "
            "two readers were NOT COMPARED. That is not the same as their "
            "agreeing; set UMAT_OTI_REPO to compare them.")
    checkout = repo / "src" / "umat_oti"
    already = sys.modules.get("umat_oti")
    if already is not None:
        imported = getattr(already, "__file__", None)
        if not runs_the_checkouts_code(imported, repo):
            where = Path(imported).resolve().parent if imported else "a namespace package"
            pytest.skip(
                f"umat_oti was already imported from {where} before this "
                f"module's tests ran, and its Python files are not "
                f"byte-identical to {checkout}, so the two readers were NOT "
                f"COMPARED against {repo}. Run this module on its own, or set "
                f"UMAT_OTI_REPO to the checkout that package was installed from.")
        return repo
    source = repo / "src"
    if str(source) not in sys.path:
        sys.path.insert(0, str(source))
    import umat_oti                                      # noqa: F401
    imported = getattr(umat_oti, "__file__", None)
    if not runs_the_checkouts_code(imported, repo):
        pytest.skip(
            f"umat_oti resolved to {imported}, which is neither {checkout} nor "
            f"a byte-identical copy of it: an install of another checkout won. "
            f"The two readers were NOT COMPARED.")
    return repo


def _store(repo: Path):
    for base in (repo.parent, repo.parent.parent):
        candidate = base / STORE_RELATIVE
        if candidate.is_file():
            return [json.loads(line) for line in
                    candidate.read_text().splitlines() if line]
    pytest.skip(f"the frozen store ({STORE_RELATIVE}) is not on this machine, "
                f"so the two readers were NOT COMPARED on real records")


@pytest.fixture(scope="module")
def both():
    repo = _producer()
    from umat_oti.contract import adapt_store_records
    rows = _store(repo)
    return rows, adapt_store_records(rows)


def test_the_two_repositories_speak_the_same_version(both):
    from umat_oti.contract import CONTRACT_VERSION as theirs
    from contract_reader import CONTRACT_VERSION as ours

    assert require_compatible(theirs, speaker="UMAT_source_transformation",
                              ours=ours) == ""
    assert theirs == ours


def test_the_two_readers_name_the_same_six_gates(both):
    from umat_oti.contract import GATES as their_gates, SEVENTH as their_seventh

    assert tuple(their_gates) == GATES
    assert their_seventh == SEVENTH


def test_the_store_row_is_not_the_contract_and_the_producer_adapts_it(both):
    """Which side does the legacy work, stated as a test.

    The frozen store predates the seventh field: no row carries
    ``primal_difference_explained_by_a_measured_control``. The producer
    recovers it from where that run did record the fact --
    ``primal.explained_by_declared_precision`` and
    ``.explained_by_operation_order`` -- and puts it in the canonical record.

    This repository does NOT do that, and must not: asking the consumer to
    re-implement the producer's legacy adapter is exactly the duplication the
    contract exists to remove, and a second implementation of a recovery rule
    is a second recovery rule. What crosses the boundary is the canonical
    record, not the row it was built from.
    """
    rows, adapted = both
    assert not any(SEVENTH in json.dumps(row) for row in rows)
    materialised = sum(
        1 for record in adapted.records
        if getattr(record.evidence, SEVENTH).is_measured())
    # 83 rows carry a measured-false primal: 19 with an explanation and 64
    # without. How many reach a RECORD depends on whether the producing
    # checkout can translate the three arguments_diverged rows, which are
    # among the 64 -- so this is 80 before those words landed and 83 after.
    assert materialised == 83 - len(adapted.errors)


def test_the_two_readers_agree_about_every_gate_in_every_record(both):
    """234 canonical records, seven fields each: 1638 three-state answers that
    must match, including which ones are NOT ESTABLISHED."""
    _, adapted = both
    from umat_oti.contract import read_gates as their_read

    compared = 0
    for record in adapted.records:
        payload = record.as_dict()
        mine = gates_of(payload)
        theirs = their_read(payload)
        for name in GATES + (SEVENTH,):
            assert mine[name].state is getattr(theirs, name).state, (
                f"{record.identity.path}: {name} is {mine[name].spelling()} "
                f"here and {getattr(theirs, name).spelling()} in the producer")
            compared += 1
    assert compared == len(adapted.records) * 7
    assert len(adapted.records) + len(adapted.errors) == 237


def test_the_two_readers_agree_about_the_settled_primal(both):
    _, adapted = both
    from umat_oti.contract import read_gates as their_read

    for record in adapted.records:
        payload = record.as_dict()
        assert primal_settled(gates_of(payload)).state is \
            their_read(payload).primal_settled().state, record.identity.path


def test_the_two_readers_agree_about_every_identity(both):
    """On the canonical record AND on the raw row, because identity is the one
    thing spelled differently in each shape and a reader that got it right in
    only one of them would silently key on the wrong field."""
    rows, adapted = both
    from umat_oti.contract import UmatIdentity

    for record in adapted.records:
        payload = record.as_dict()
        assert identity_of(payload) == (record.identity.path,
                                        record.identity.sha256)
    for row in rows:
        theirs = UmatIdentity.from_record(row)
        assert identity_of(row) == (theirs.path, theirs.sha256)


def test_the_two_readers_agree_about_every_terminal_state_and_owner(both):
    rows, adapted = both
    by_path = {r.identity.path: r for r in adapted.records}
    refused = {e["identity"]["path"] for e in adapted.errors}

    for row in rows:
        path, _ = identity_of(row)
        if path in refused:
            # The producer refused to translate the stage. This reader must
            # refuse the resulting state too, rather than inventing an owner.
            assert path not in by_path
            continue
        theirs = by_path[path]
        assert owner_of(theirs.terminal.state) == theirs.terminal.owner


def test_a_stage_nobody_mapped_is_refused_and_never_becomes_not_attempted(both):
    """Three of the 237 entries carry a rung whose word arrived late.

    They RAN, and the vocabulary's old fall-through answered
    ``not_attempted`` -- INTERNAL, "the run never happened" -- about all
    three, booking an external fact into the producing project's column.

    Two states of the world are acceptable and this accepts either, because
    which one holds depends on whether the producing checkout has the word
    yet. What is never acceptable is ``not_attempted``.
    """
    rows, adapted = both
    diverged = [r for r in rows
                if r["stage"] == "arguments_diverged_before_the_routine"]
    assert len(diverged) == 3

    if adapted.errors:
        # The producer's FROM_STAGE has no entry, so it refuses rather than
        # defaulting, and every refusal is accounted for.
        assert len(adapted.errors) == 3
        for error in adapted.errors:
            assert error["code"] == "contract.untranslatable_stage"
            assert error["detail"]["stage"] == \
                "arguments_diverged_before_the_routine"
            assert "not_attempted" in error["message"]
    else:
        # The producer has the word. It is the producer's own (INTERNAL), and
        # this repository must know it too or it would fall through.
        states = {r.terminal.state for r in adapted.records
                  if r.terminal.stage == "arguments_diverged_before_the_routine"}
        assert states == {"arguments_diverged_before_the_routine"}
        assert owner_of("arguments_diverged_before_the_routine") == "INTERNAL"

    # Either way, nothing anywhere became not_attempted.
    assert all(r.terminal.state != "not_attempted" for r in adapted.records)
    # And a rung nothing has ever mapped is still refused on this side.
    with pytest.raises(ContractError):
        owner_of("a_rung_nobody_has_ever_mapped")


def test_the_two_readers_agree_about_the_whole_terminal_vocabulary(both):
    """The consuming repository cannot import the producer's vocabulary, so it
    carries the published one. Any state the producer can emit that this
    reader cannot own would arrive here and fall through."""
    from umat_oti.contract.terminal import PUBLISHED_OWNERS

    for state, owner in PUBLISHED_OWNERS.items():
        assert owner_of(state) == owner, state


def test_the_two_readers_agree_about_which_entries_may_drive_the_assembler(both):
    """The answer the assembler actually acts on, on all 234 translated
    records -- and three-state, so "nobody measured this" cannot come out as
    "go ahead" on either side."""
    _, adapted = both
    usable = 0
    for record in adapted.records:
        payload = record.as_dict()
        mine, why = usable_for_assembly(payload)
        theirs, their_why = record.usable_by_the_residual_assembler()
        assert mine.state is theirs.state, (
            f"{record.identity.path}: {mine.spelling()} here, "
            f"{theirs.spelling()} in the producer ({why!r} / {their_why!r})")
        if mine.is_true():
            usable += 1
            assert record.terminal.verified
    # 55 is fixed: it is how many entries reached a verified terminal state,
    # and no vocabulary change moves it. The denominator is what moves.
    assert usable == 55
    assert len(adapted.records) + len(adapted.errors) == 237


def test_no_record_is_usable_without_all_six_gates_settled(both):
    """The one-line summary of the whole contract, checked on real data."""
    _, adapted = both
    for record in adapted.records:
        answer, _ = usable_for_assembly(record.as_dict())
        if not answer.is_true():
            continue
        gates = gates_of(record.as_dict())
        assert primal_settled(gates).is_true()
        for name in GATES:
            if name == "primal_agreed":
                continue
            assert gates[name].is_true(), f"{record.identity.path}: {name}"
