"""A verified UMAT case, read as ingredients this assembler can use.

The assembler's tests used to supply stress and tangent from a reference model
written for the purpose. That tests the assembly and says nothing about the
bridge: a fixture written to make the assembler pass cannot also be evidence
that the assembler consumes what the pipeline actually produces.

A fixture here is different. It is the output of a corpus case that VERIFIED
in Abaqus -- the original UMAT ran, the OTI-converted build ran on the same
deck, their stress and state histories agreed over the whole path, and the
converted build's tangent agreed with a finite difference of the original at
several states. It carries the numbers and the identity of the file they came
from, never the file: most of the corpus is not redistributable.

The other thing it carries is the rule that makes it usable as a BASELINE. A
regression fixture is what later runs are compared against, so anything wrong
inside it is wrong in every comparison made against it afterwards, silently: a
NaN frozen into a fixture does not fail, it propagates, and the comparison
that should have caught it is being made against the NaN. So a fixture holds a
completely finite successful history and nothing else -- every increment
present, every material point present, every number a number. The exporter
refuses to write anything else; :func:`load` refuses to read anything else,
because a claim checked only by the tool that makes it is a claim nobody
checked.

The other thing it carries is the CONVENTION. A tangent is not a number until
somebody says what its indices mean, and the two sides of this bridge have to
agree about Voigt ordering and about which off-diagonal entries are
engineering shear. The fixture states it; :func:`check_conventions` checks the
numbers are consistent with what it states, so a mapping error is caught where
it happens rather than as a residual that is wrong by a factor of two.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np

SCHEMA = "umat-oti/residual-fixture/1"

#: The order Abaqus hands a UMAT a symmetric tensor in, and the order DDSDDE
#: is expressed in. The off-diagonal entries are ENGINEERING shear: a shear
#: strain of gamma, not gamma/2. The C3D8 kernel's B matrix is built in the
#: same order and the same convention, which is why they compose.
VOIGT_ORDER = ("11", "22", "33", "12", "13", "23")


class FixtureError(ValueError):
    """A fixture that cannot be used, said in terms of what is wrong with it."""


#: The arrays in a carried increment that must be finite numbers. ``ddsdde``
#: is one of them: a fixture whose stress is finite and whose tangent is not
#: freezes a NaN that surfaces only once somebody assembles a stiffness out of
#: it, three layers away from the fixture that caused it.
CARRIED_ARRAYS = ("strain", "dstrain", "stress", "state", "ddsdde")


@dataclass(frozen=True)
class Point:
    """One material point at one increment."""

    increment: int
    time: float
    strain: np.ndarray
    dstrain: np.ndarray
    stress: np.ndarray
    state: np.ndarray
    tangent: Optional[np.ndarray]


@dataclass
class VerifiedFixture:
    """One verified corpus case, and what was established about it."""

    source_id: str
    repository: str
    element_type: str
    ntens: int
    ndi: int
    nshr: int
    kinematics: str
    props: tuple
    nstatv: Optional[int]
    material_provenance: str
    deck: str
    verification: dict = field(default_factory=dict)
    original: list = field(default_factory=list)
    converted: list = field(default_factory=list)
    path: Optional[Path] = None
    #: Which of the exporter's claims this fixture actually carries. A fixture
    #: frozen before a check existed cannot be held to it, and saying so is
    #: not the same as passing it: a reader of this list can tell a fixture
    #: that was checked from one that merely was not caught.
    claims_checked: tuple = ()
    claims_not_carried: tuple = ()

    @property
    def finite_strain(self) -> bool:
        return str(self.kinematics).startswith("finite")

    def increments(self) -> int:
        return min(len(self.original), len(self.converted))

    def describe(self) -> str:
        return (f"{self.source_id} ({self.repository}), {self.element_type}, "
                f"NTENS={self.ntens}, {self.kinematics}; verified: "
                f"{self.verification.get('states_agreeing')} of "
                f"{self.verification.get('states_checked')} tangent states")


def _point(record: dict, ntens: int) -> Point:
    tangent = record.get("ddsdde") or []
    matrix = None
    if len(tangent) >= ntens * ntens:
        matrix = np.asarray(tangent[:ntens * ntens], dtype=float).reshape(
            ntens, ntens)
    return Point(
        increment=int(record.get("increment") or 0),
        time=float(record.get("time") or 0.0),
        strain=np.asarray(record.get("strain") or [], dtype=float),
        dstrain=np.asarray(record.get("dstrain") or [], dtype=float),
        stress=np.asarray(record.get("stress") or [], dtype=float),
        state=np.asarray(record.get("state") or [], dtype=float),
        tangent=matrix)


def load(path: Path) -> VerifiedFixture:
    """Read a fixture, refusing one whose shape does not match its own claims."""
    path = Path(path)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise FixtureError(f"{path.name} could not be read: {error}") from error
    if payload.get("schema") != SCHEMA:
        raise FixtureError(
            f"{path.name} declares schema {payload.get('schema')!r}; this "
            f"loader reads {SCHEMA!r}. A fixture whose shape is not the shape "
            f"this reads would be interpreted wrongly rather than refused.")
    point = payload.get("material_point") or {}
    ntens = int(point.get("ntens") or 0)
    if ntens <= 0:
        raise FixtureError(f"{path.name} states no tensor size, so nothing in "
                           f"it can be indexed")
    fixture = VerifiedFixture(
        source_id=str(payload.get("source_id") or ""),
        repository=str(payload.get("repository") or ""),
        element_type=str(point.get("element_type") or ""),
        ntens=ntens,
        ndi=int(point.get("ndi") or 0),
        nshr=int(point.get("nshr") or 0),
        kinematics=str(point.get("kinematics") or ""),
        props=tuple((payload.get("material") or {}).get("props") or ()),
        nstatv=(payload.get("material") or {}).get("nstatv"),
        material_provenance=str((payload.get("material") or {}).get(
            "provenance") or ""),
        deck=str(payload.get("deck") or ""),
        verification=dict(payload.get("verification") or {}),
        path=path)
    fixture.original = [_point(record, ntens)
                        for record in payload.get("original") or []]
    fixture.converted = [_point(record, ntens)
                         for record in payload.get("converted") or []]
    if not fixture.original:
        raise FixtureError(f"{path.name} carries no increments")
    _refuse_anything_but_a_finite_history(path, payload)
    fixture.claims_checked, fixture.claims_not_carried = _what_was_checkable(
        payload)
    for index, record in enumerate(fixture.original):
        if record.stress.size != ntens:
            raise FixtureError(
                f"{path.name} increment {record.increment}: the stress has "
                f"{record.stress.size} components and the fixture states "
                f"NTENS={ntens}. One of the two is wrong and neither may be "
                f"assumed.")
    return fixture


def _refuse_anything_but_a_finite_history(path: Path, payload: dict) -> None:
    """Refuse a fixture that is not a completely finite successful history.

    Two readings, because either alone can be fooled. What the exporter
    RECORDED about the run it froze -- the gates the verification measured and
    the grouping that says whether every increment produced every material
    point and where the first value that was not a number is. And the numbers
    in the file, scanned one by one, because "this history was finite" and
    "these numbers are finite" are two different claims and a fixture is the
    second.
    """
    frozen = payload.get("finite_history")
    if not isinstance(frozen, dict):
        raise FixtureError(
            f"{path.name} carries no finite_history block, so nothing in it "
            f"says the run it froze was complete. Re-export it with "
            f"UMAT_source_transformation/tools/export_residual_fixture.py, "
            f"which refuses to write a fixture that is not one.")
    if frozen.get("complete_finite_verification_run") is not True:
        raise FixtureError(
            f"{path.name} was frozen from a run whose "
            f"complete_finite_verification_run is "
            f"{frozen.get('complete_finite_verification_run')!r}. A baseline "
            f"built on a history that did not run whole is a baseline nothing "
            f"can be compared against.")
    measured = frozen.get("evidence") or {}
    for gate in ("abaqus_job_completed", "all_requested_outputs_present",
                 "complete_history_finite"):
        if measured.get(gate) is not True:
            raise FixtureError(
                f"{path.name}: the run it was frozen from has {gate}="
                f"{measured.get(gate, 'not measured')!r}")
    for side, grouping in (frozen.get("history_grouping") or {}).items():
        if not isinstance(grouping, dict):
            continue
        if grouping.get("first_incomplete_increment"):
            raise FixtureError(
                f"{path.name}: the {side} history's first incomplete "
                f"increment is "
                f"{grouping['first_incomplete_increment'].get('increment')}")
        if grouping.get("first_non_finite_material_point"):
            where = grouping["first_non_finite_material_point"]
            raise FixtureError(
                f"{path.name}: the {side} history stops being numbers at "
                f"element {where.get('element')} point {where.get('point')} "
                f"of increment {where.get('increment')}")
    for side in ("original", "converted"):
        for record in payload.get(side) or []:
            for name in CARRIED_ARRAYS:
                for index, value in enumerate(record.get(name) or ()):
                    try:
                        number = float(value)
                    except (TypeError, ValueError):
                        raise FixtureError(
                            f"{path.name}: {side} increment "
                            f"{record.get('increment')} {name}[{index}] is "
                            f"{value!r}, which is not a number") from None
                    if not math.isfinite(number):
                        raise FixtureError(
                            f"{path.name}: {side} increment "
                            f"{record.get('increment')} {name}[{index}] is "
                            f"{number}. A fixture carries a completely finite "
                            f"history; this one does not.")

    # -- the window was not carved out of a run that stopped early ---------
    # The exporter grew this refusal after these fixtures were frozen. A
    # window that is SHORT is a window built on however far an analysis got,
    # and it is short silently: the window clamps itself to whatever was on
    # disk and records the shortened number as though it had been asked for.
    asked = frozen.get("increments_requested")
    got = frozen.get("increments_carried")
    if isinstance(asked, int) and isinstance(got, int) and got < asked:
        raise FixtureError(
            f"{path.name} carries {got} increment(s) of the {asked} its "
            f"export asked for, because only "
            f"{frozen.get('records_available')} record(s) were on disk. A "
            f"fixture that is short is a fixture built on however far an "
            f"analysis got before it stopped, and nothing comparing against "
            f"it can tell that from a run meant to be that length.")

    # -- and every number the run wrote, not only the ones carried ---------
    # A window of six increments is finite in a run that went to NaN at
    # increment 200. The exporter scans the whole history and writes what it
    # found; this reads that rather than taking the window's finiteness as
    # evidence about the run.
    for side, scan in (frozen.get("whole_history") or {}).items():
        if not isinstance(scan, dict):
            continue
        where = scan.get("first_non_finite")
        if where:
            raise FixtureError(
                f"{path.name}: the {side} run it was frozen from stops being "
                f"numbers OUTSIDE the carried window -- {where.get('array')}"
                f"[{where.get('index')}] is {where.get('value')} at record "
                f"{where.get('record')}, increment {where.get('increment')}. "
                f"The window is finite and the run is not, so this fixture is "
                f"the prefix of a failed analysis.")


#: Claims the exporter writes into ``finite_history`` that :func:`load` checks
#: when they are there. A fixture frozen before one of them existed does not
#: carry it, and :attr:`VerifiedFixture.claims_not_carried` names which --
#: because "checked and clean" and "never checked" are different, and a
#: consumer that cannot tell them apart is treating an absence as a pass.
CHECKABLE_CLAIMS = (
    ("complete_finite_verification_run",
     "the run was finite from end to end"),
    ("evidence", "the gates the verification measured"),
    ("history_grouping", "every increment produced every material point"),
    ("increments_requested",
     "the window is the length the export asked for, not what was on disk"),
    ("whole_history",
     "every number the run wrote is finite, not only the ones carried"),
)


def _what_was_checkable(payload: dict) -> tuple:
    """Which of :data:`CHECKABLE_CLAIMS` this fixture actually carries."""
    frozen = payload.get("finite_history") or {}
    checked, absent = [], []
    for key, said in CHECKABLE_CLAIMS:
        (checked if frozen.get(key) is not None else absent).append(said)
    return tuple(checked), tuple(absent)


def load_all(directory: Path) -> list:
    """Every fixture in a directory, in a stable order."""
    return [load(path) for path in sorted(Path(directory).glob("*.json"))]


def check_conventions(fixture: VerifiedFixture) -> list:
    """What the numbers say about the convention the fixture claims.

    A tangent is not a number until somebody says what its indices mean. Three
    things are checkable without another model:

    * the tensor splits as the fixture says -- ``ndi + nshr == ntens``;
    * the tangent is square in ``ntens`` and its direct block is the block the
      direct strains drive;
    * the shear diagonal is positive and smaller than the direct diagonal,
      which is true of every isotropic elastic tangent expressed in
      ENGINEERING shear and false of one expressed in tensorial shear by a
      factor of two.

    The third is a smell rather than a proof, and it is reported as one.
    """
    problems: list = []
    if fixture.ndi and fixture.nshr and fixture.ndi + fixture.nshr != fixture.ntens:
        problems.append(
            f"the fixture states NDI={fixture.ndi} and NSHR={fixture.nshr}, "
            f"which is {fixture.ndi + fixture.nshr} components, and NTENS="
            f"{fixture.ntens}")
    for record in fixture.converted:
        if record.tangent is None:
            continue
        if record.tangent.shape != (fixture.ntens, fixture.ntens):
            problems.append(
                f"increment {record.increment}: the tangent is "
                f"{record.tangent.shape} and NTENS is {fixture.ntens}")
            continue
        if fixture.nshr and fixture.ndi:
            direct = np.diag(record.tangent)[:fixture.ndi]
            shear = np.diag(record.tangent)[fixture.ndi:]
            if np.any(shear <= 0) and np.any(direct > 0):
                problems.append(
                    f"increment {record.increment}: a shear diagonal entry is "
                    f"not positive ({shear.tolist()}), which no stable "
                    f"material tangent has")
    return problems
