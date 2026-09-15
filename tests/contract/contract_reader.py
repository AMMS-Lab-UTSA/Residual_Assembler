"""What this repository needs in order to READ a UMAT contract record.

Deliberately small, and deliberately free of ``residual_core``. It is the
consumer half of the agreement: version handshake, three-state gate reading,
identity, terminal-state ownership, and the fixture fingerprint check. The
producing repository has a fuller implementation; this one exists so that
Residual_Assembler can state what it speaks and refuse what it cannot read
without depending on the other repository being installed.

The three-state reader itself is NOT reimplemented here. It is imported from
``_shared/tristate.py``, which is a byte-identical copy of the producing
repository's module, held to that by ``schemas/contract_lock.json``. A
reimplementation is exactly where "null reads as a pass" comes back.
"""
from __future__ import annotations

import json
import posixpath
import re
from pathlib import Path
from typing import Any, Mapping

from frames import (COUNT_FIELDS, FrameError, IDENTITY_FIELDS,  # noqa: F401
                    check_counts, count_history, frame_key,
                    group_by_increment, increment_key, point_key)
from tristate import NOT_ESTABLISHED, Tri, all_true, read  # noqa: F401

#: Which shared data contract THIS repository speaks. A mismatch with the
#: producer is refused loudly rather than surfacing as a missing key.
#:
#: 2.0.0: two terminal states this repository had never heard of, a row
#: identity that grew from two fields to five, five counts where there was
#: one, and a result envelope whose success is no longer shaped like its
#: verdict. A 1.x reader would have MISREAD every one of those rather than
#: failed on them, which is what makes the bump major.
CONTRACT_VERSION = "2.0.0"
SPEAKER = "Residual_Assembler"

SHA256 = re.compile(r"^[0-9a-f]{64}$")

#: The six gates, in the order the contract reports them. Six, never one.
GATES = (
    "abaqus_job_completed",
    "all_requested_outputs_present",
    "complete_history_finite",
    "primal_agreed",
    "derivatives_verified",
    "mechanically_informative",
)
#: The seventh field. Not a gate: it says why a false primal is still a verdict.
SEVENTH = "primal_difference_explained_by_a_measured_control"

#: A fact about somebody's published repository. More engineering in the
#: producing project changes none of them.
#:
#: ``arguments_diverged_before_the_routine`` is here and is EXTERNAL: the
#: paired call's INPUTS already differed before the routine was entered, so
#: the difference is upstream of the constitutive code. Three of the 237
#: frozen entries carry it, and before the vocabulary had a word for them they
#: were answered ``not_attempted`` -- INTERNAL, "the run never happened" --
#: about three entries that ran.
EXTERNAL_STATES = frozenset({
    "missing_material_data", "not_a_umat", "incomplete_or_corrupt_source",
    "external_dependency_unavailable",
    "published_stub_no_constitutive_content", "waits_for_input",
    "arguments_diverged_before_the_routine"})
#: A limitation of the producing project. Every one of these is work there.
#:
#: ``primal_mismatch_explained`` is here and is INTERNAL: the primal results
#: disagreed and a measured control accounts for the difference. An explained
#: disagreement is still a disagreement -- the primal gate is false -- so it is
#: neither verified nor a fact about somebody else's repository.
INTERNAL_STATES = frozenset({
    "transform_refused", "experiment_not_generated",
    "experiment_not_informative", "informativeness_not_established",
    "unsupported_formulation", "support_build_failed", "original_job_failed",
    "transformed_job_failed", "primal_disagreed", "primal_mismatch_explained",
    "disagreement_not_in_any_recorded_call", "derivative_truncated",
    "tangent_not_verified", "not_attempted", "harness_error"})
FULLY_VERIFIED = "fully_verified"


class ContractError(RuntimeError):
    """This repository will not read what it has been handed."""


class ContractVersionError(ContractError):
    """The two ends do not speak compatible versions of the contract."""


def require_compatible(theirs: Any, *, speaker: str = "the producer",
                       ours: str = CONTRACT_VERSION) -> str:
    """Refuse an incompatible producer at the boundary, with both versions named."""
    def parse(text: Any) -> tuple:
        if text is None:
            raise ContractVersionError(
                f"this document declares no contract version. {SPEAKER} "
                f"speaks {ours} and cannot assume the fields mean what {ours} "
                f"says they mean.")
        parts = str(text).split(".")
        if len(parts) != 3 or not all(p.isdigit() for p in parts):
            raise ContractVersionError(
                f"{text!r} is not a contract version; expected "
                f"'MAJOR.MINOR.PATCH', for example {ours!r}.")
        return tuple(int(p) for p in parts)

    mine, yours = parse(ours), parse(theirs)
    if yours[0] != mine[0]:
        raise ContractVersionError(
            f"contract version mismatch: {speaker} speaks "
            f"{'.'.join(map(str, yours))}, {SPEAKER} speaks "
            f"{'.'.join(map(str, mine))}. A major-version difference is a "
            f"BREAKING change -- fields removed, renamed, retyped or given a "
            f"new meaning -- so reading this would not fail, it would produce "
            f"wrong numbers. Refusing.")
    if yours[1] > mine[1]:
        raise ContractVersionError(
            f"contract version mismatch: {speaker} speaks "
            f"{'.'.join(map(str, yours))}, newer than {SPEAKER}'s "
            f"{'.'.join(map(str, mine))}. This reader cannot tell a field it "
            f"has never heard of from a field that is absent. Refusing rather "
            f"than reading it partially.")
    if yours[1] < mine[1]:
        return (f"{speaker} speaks {'.'.join(map(str, yours))}, older than "
                f"{SPEAKER}'s {'.'.join(map(str, mine))}; every field added "
                f"since is NOT ESTABLISHED in this record.")
    return ""


def handshake() -> dict:
    return {"speaker": SPEAKER, "contract_version": CONTRACT_VERSION}


# ---------------------------------------------------------------------------
def identity_of(record: Mapping[str, Any]) -> tuple:
    """``(path, sha256)``. Refuses a bare basename: 23 basenames in the frozen
    corpus name more than one source file, covering 108 of 237 entries.

    Three shapes carry an identity and this reads all three, because the
    producer's reader does and two readers that disagree about where the
    identity lives are two contracts. A contract record puts it in an
    ``identity`` block; a frozen fixture spells it ``source_id``; a raw
    store-verification row spells it ``source``. The DIGEST field is spelled
    ``source_sha256`` in both of the flat shapes.
    """
    block = record.get("identity") if isinstance(record, Mapping) else None
    if isinstance(block, Mapping):
        path, digest = block.get("path"), block.get("sha256")
    else:
        path = record.get("source_id") or record.get("source")
        digest = record.get("source_sha256")
    text = str(path or "").replace("\\", "/").strip()
    if not text:
        raise ContractError("this record carries no source path, so it cannot "
                            "be matched to a source at all")
    if text.startswith("/"):
        raise ContractError(f"{text!r} is absolute; identity is the path "
                            f"WITHIN the corpus cache")
    if posixpath.dirname(text) == "":
        raise ContractError(
            f"{text!r} is a bare basename, which is not an identity. A "
            f"consumer keyed on one silently merges materials that have "
            f"different constants and different verdicts.")
    digest = str(digest or "").strip().lower()
    if not SHA256.match(digest):
        raise ContractError(
            f"{digest!r} is not a SHA-256. A path without a digest says which "
            f"file was meant, not which bytes were verified.")
    return (text, digest)


def gates_of(record: Mapping[str, Any]) -> dict:
    """The seven, three-state. A missing block leaves every one NOT ESTABLISHED."""
    evidence = record.get("evidence")
    if evidence is not None and not isinstance(evidence, Mapping):
        raise ContractError(
            f"'evidence' must be an object of gate fields or absent; found "
            f"{type(evidence).__name__}, which cannot be read three-state")
    return {name: read(evidence, name) for name in GATES + (SEVENTH,)}


def primal_settled(gates: Mapping[str, Tri]) -> Tri:
    """The primal after its control -- the chain a verdict rests on."""
    primal, seventh = gates["primal_agreed"], gates[SEVENTH]
    if primal.is_true():
        return Tri(True, "the two builds agreed")
    if primal.is_not_established():
        return Tri(None, "nothing compared the two builds")
    if seventh.is_true():
        return Tri(True, "a measured control explained the difference")
    if seventh.is_false():
        return Tri(False, "a control ran and did not explain the difference")
    return Tri(None, "the builds disagreed and nothing says whether a control "
                     "explained it; an unexplained disagreement is not a "
                     "verification")


def owner_of(state: str) -> str:
    """EXTERNAL, INTERNAL or NONE. An unknown state is refused, not defaulted."""
    if state == FULLY_VERIFIED:
        return "NONE"
    if state in EXTERNAL_STATES:
        return "EXTERNAL"
    if state in INTERNAL_STATES:
        return "INTERNAL"
    raise ContractError(
        f"{state!r} is not a terminal state {SPEAKER} knows. It cannot say "
        f"whether this is EXTERNAL -- a fact about somebody's published "
        f"repository -- or INTERNAL, a limitation of the producing project, "
        f"and guessing would assert an owner nobody established.")


def check_fingerprint(fixture_fingerprint: Any, store_fingerprint: Any, *,
                      where: str = "this fixture") -> Tri:
    """Three-state: current, stale, or carrying no fingerprint at all."""
    theirs = str(fixture_fingerprint or "").strip()
    mine = str(store_fingerprint or "").strip()
    if not theirs:
        return Tri(None, f"{where} carries no transform_fingerprint, so "
                         f"whether its numbers predate a correction to the "
                         f"transform is not established -- it is not current "
                         f"by default")
    if not mine:
        return Tri(None, f"no store fingerprint was supplied for {where}")
    if theirs == mine:
        return Tri(True, f"{where} was frozen at {mine}")
    return Tri(False,
               f"{where} was frozen at transform fingerprint {theirs} and the "
               f"store it is read beside was produced at {mine}. Its stress, "
               f"state and tangent were computed by a transform that has "
               f"since changed, so differencing against them is evidence "
               f"about the old transform and not the current one.")


def require_current(fixture_fingerprint: Any, store_fingerprint: Any, *,
                    where: str = "this fixture") -> None:
    """Raise on anything but TRUE. NOT ESTABLISHED is not known-current."""
    answer = check_fingerprint(fixture_fingerprint, store_fingerprint,
                               where=where)
    if not answer.is_true():
        raise ContractError(answer.why)


def usable_for_assembly(record: Mapping[str, Any]) -> tuple:
    """``(Tri, reason)`` -- may the assembler drive itself from this entry?

    Three-state. A refusal is a finding; an entry nobody finished measuring is
    a queue item; and merging them would put both into one bucket when the
    work they imply is different.
    """
    terminal = record.get("terminal") or {}
    state = str(terminal.get("state") or "")
    if state != FULLY_VERIFIED:
        owner = owner_of(state) if state else "unknown"
        return (Tri(False), f"terminal state is {state!r} ({owner}), not "
                            f"{FULLY_VERIFIED}")
    gates = gates_of(record)
    settled = all_true(primal_settled(gates),
                       *(gates[n] for n in GATES if n != "primal_agreed"))
    primal = primal_settled(gates)
    if settled.is_false():
        failing = [n for n in GATES
                   if gates[n].is_false() and n != "primal_agreed"]
        if primal.is_false():
            failing.append(f"primal_agreed ({primal.why})")
        return (Tri(False), "a gate was measured false: " + ", ".join(failing))
    if settled.is_not_established():
        unmeasured = [n for n in GATES
                      if gates[n].is_not_established() and n != "primal_agreed"]
        if primal.is_not_established():
            unmeasured.append(f"primal_agreed ({primal.why})")
        return (NOT_ESTABLISHED,
                "a gate was never measured: " + ", ".join(unmeasured))
    convention = (record.get("convention") or {}).get("voigt_order") or []
    if not convention:
        return (NOT_ESTABLISHED, "the tensor component order is not "
                                 "established, so an assembled residual "
                                 "could not be indexed")
    return (Tri(True), "")


# ---------------------------------------------------------------------------
# an envelope is not a verdict
# ---------------------------------------------------------------------------
#: Field names this repository must never read as a verdict about a material,
#: with what each actually means. ``ok`` is the one that did the damage: it
#: says the service RAN, which is true of every entry it managed to look at
#: including the ones it refuses, and reading it as the answer offered all 237
#: corpus materials here as verified when 42 are true on all six gates.
FORBIDDEN_AS_A_VERDICT = {
    "ok": "whether the producing service RAN, not what it found",
    "success": "whether the call succeeded, not what it found",
    "passed": "ambiguous between the call and the subject",
    "status": "a service's own word for what happened",
}
SUCCESS_FIELD = "call_succeeded"
VERDICT_FIELD = "verdict"


def verdict_of(payload: Mapping[str, Any]) -> Tri:
    """The verdict a payload states, refusing to accept a call-success for one.

    Raises rather than shrugging when the payload states none: answering NOT
    ESTABLISHED would hide the misreading, and the caller asked a question
    this payload does not answer.
    """
    if not isinstance(payload, Mapping):
        raise ContractError(f"expected an object; got {type(payload).__name__}")
    if VERDICT_FIELD in payload:
        value = payload[VERDICT_FIELD]
        if value is None or isinstance(value, bool):
            return Tri(value)
        raise ContractError(
            f"{VERDICT_FIELD} must be true, false or null; got {value!r}")
    present = [n for n in FORBIDDEN_AS_A_VERDICT if n in payload]
    if present:
        raise ContractError(
            f"this payload states no {VERDICT_FIELD!r}; it carries "
            f"{', '.join(repr(n) for n in present)} and none of those is one ("
            + "; ".join(f"{n}: {FORBIDDEN_AS_A_VERDICT[n]}" for n in present)
            + "). Reading one of them as the answer offered 237 materials as "
              "verified when 42 are true on all six gates.")
    raise ContractError(
        f"this payload states no {VERDICT_FIELD!r} and carries nothing that "
        f"could be mistaken for one. Supply the verdict explicitly.")


# ---------------------------------------------------------------------------
# which generation a fixture speaks
# ---------------------------------------------------------------------------
def fixture_generation(payload: Mapping[str, Any]) -> dict:
    """Whether a fixture's rows can be told apart, and its counts trusted.

    A fixture carries no contract version -- it predates the contract -- so
    the generation is read off its shape. A 1.x fixture is not broken and is
    not current either: its numbers may be differenced against, and the
    loading behind them may not be reconstructed from it.
    """
    rows = payload.get("original") or []
    first = rows[0] if rows and isinstance(rows[0], Mapping) else {}
    five = all(f in first for f in IDENTITY_FIELDS)
    history = payload.get("finite_history")
    history = history if isinstance(history, Mapping) else {}
    counts = all(history.get(n) is not None for n in
                 ("records_carried", "increments_carried",
                  "material_points_per_increment"))
    absent = [f for f in IDENTITY_FIELDS if f not in first]
    return {
        "identity": "2.x" if five else "1.x",
        "has_counts": counts,
        "usable_for_boundary_conditions": five,
        "reason": ("rows carry all five identity fields" if five else
                   f"rows carry no {', '.join(absent)}; Abaqus restarts "
                   f"increment numbering in every step, so the loading cannot "
                   f"be rebuilt from what is here -- doing it anyway was wrong "
                   f"by 3.0 relative, a different deformation and not a "
                   f"tolerance. The numbers may still be differenced against."),
    }


def require_five_field_identity(payload: Mapping[str, Any], *,
                                where: str = "this fixture") -> None:
    """Refuse a fixture whose rows cannot be told apart, for a consumer that
    has to rebuild the loading. A consumer that only differences stored
    numbers must not call this: refusing a usable fixture is its own wrong
    answer."""
    generation = fixture_generation(payload)
    if not generation["usable_for_boundary_conditions"]:
        raise ContractError(f"{where}: {generation['reason']}")
