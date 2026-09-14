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

#: The transform fingerprint of the store these fixtures are allowed to have
#: come from. A fixture is a number the OTI transformation produced, and the
#: transformation changes: the same UMAT converted by a different build of it
#: is a different experiment, and its stress and tangent are evidence about
#: that build rather than about this one. Verifying today's assembler against
#: a fixture frozen under an older fingerprint is verifying it against
#: somebody else's run.
#:
#: Measured, not assumed: every one of the 237 entries in the store named
#: below carries this fingerprint, and the four fixtures this repository
#: carried before 2026-09-14 carried ``ff94800b1884bcc0`` -- a store two
#: transformations ago.
#: READ FROM THE CONTRACT, not written here.
#:
#: It was a literal, and it went stale the moment the transformation changed:
#: the store moved to a new build, every fixture was re-frozen against it, and
#: this reader went on accepting the old value and refusing all ten of them.
#: A number that two repositories must agree on belongs in the one file they
#: both read -- ``schemas/transform_generation.json``, which the UMAT
#: repository writes and copies here, and whose own ``how_to_update`` describes
#: this exact sequence.
def _recorded_generation() -> str:
    """The transform the current evidence was cut at."""
    import json as _json

    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "schemas" / "transform_generation.json"
        if candidate.is_file():
            recorded = _json.loads(candidate.read_text(encoding="utf-8"))
            value = str(recorded.get("transform_fingerprint") or "")
            if value:
                return value
    raise FixtureError(
        "schemas/transform_generation.json is missing, so nothing here knows "
        "which build of the transformation the committed evidence was cut at. "
        "It is copied from the UMAT repository; without it a fixture cannot be "
        "told from evidence about a transformation that no longer exists.")


CURRENT_TRANSFORM_FINGERPRINT = _recorded_generation()

#: Where that fingerprint was read from, so a reader can check it rather than
#: take it. Counts are from the file itself.
STORE_PROVENANCE = {
    "results": ("corpus_run/pass12/results/store_verification.jsonl"),
    "entries": 244,
    "verified": 57,
    "all_six_evidence_gates": 44,
    "fingerprint": CURRENT_TRANSFORM_FINGERPRINT,
    "read_from": "schemas/transform_generation.json",
    "read_on": "2026-09-14",
}

#: The six gates the UMAT pipeline measures on a case. A fixture carries what
#: was measured; it is the reader's business what to require. ``derivatives_
#: verified`` and ``primal_agreed`` are the two that can be false on a case
#: that still reached stage "verified", and a consumer that does not look at
#: them is treating "verified" as a single bit when the pipeline recorded six.
EVIDENCE_GATES = ("abaqus_job_completed", "all_requested_outputs_present",
                  "complete_history_finite", "derivatives_verified",
                  "mechanically_informative", "primal_agreed")

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
    """One material point at one increment.

    ``step`` and ``point`` complete the identity. Abaqus numbers increments
    from 1 again in every step, so ``increment`` alone does not name a state:
    a four-step cycle has four increment 1s, and reconstructing the strain
    without the step compared a reversal against the first loading and was
    wrong by 3.0 relative. Both are None in a fixture frozen before the
    exporter carried them, and those are the single-step, single-point ones.
    """

    increment: int
    time: float
    strain: np.ndarray
    dstrain: np.ndarray
    stress: np.ndarray
    state: np.ndarray
    tangent: Optional[np.ndarray]
    step: Optional[int] = None
    point: Optional[int] = None
    element: Optional[int] = None


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
    #: The transform build that produced these numbers, and the six gates the
    #: pipeline measured on the run they came from.
    transform_fingerprint: str = ""
    source_sha256: str = ""
    deck_digest: str = ""
    evidence: dict = field(default_factory=dict)

    @property
    def finite_strain(self) -> bool:
        return str(self.kinematics).startswith("finite")

    @property
    def gates_true(self) -> tuple:
        return tuple(g for g in EVIDENCE_GATES if self.evidence.get(g) is True)

    @property
    def gates_not_true(self) -> tuple:
        """Gates that are false OR were never measured, kept apart from the
        ones that passed. A gate absent from the record is not a gate that
        held."""
        return tuple(g for g in EVIDENCE_GATES if self.evidence.get(g) is not True)

    @property
    def all_six_gates(self) -> bool:
        return not self.gates_not_true

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
        tangent=matrix,
        step=(int(record["step"]) if record.get("step") is not None else None),
        point=(int(record["point"]) if record.get("point") is not None else None),
        element=(int(record["element"]) if record.get("element") is not None
                 else None))


def load(path: Path, *,
         fingerprint: Optional[str] = CURRENT_TRANSFORM_FINGERPRINT
         ) -> VerifiedFixture:
    """Read a fixture, refusing one whose shape does not match its own claims.

    ``fingerprint`` is the transform build this reader accepts evidence from;
    a fixture frozen under a different one is refused, because the numbers in
    it describe a transformation that is no longer the one being tested.
    Passing ``None`` reads a fixture of any provenance -- which is a
    deliberate act with a reason, not a default.
    """
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
    carried = str(payload.get("transform_fingerprint") or "")
    if fingerprint is not None and carried != fingerprint:
        raise FixtureError(
            f"{path.name} was frozen under transform fingerprint "
            f"{carried or 'none recorded'!r} and this reader accepts "
            f"{fingerprint!r} ({STORE_PROVENANCE['results']}). The same UMAT "
            f"converted by a different build of the transformation is a "
            f"different experiment: its stress and tangent are evidence about "
            f"that build, so checking today's assembler against them would be "
            f"checking it against somebody else's run. Re-export it from the "
            f"current store with "
            f"UMAT_source_transformation/tools/export_residual_fixture.py.")
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
        path=path,
        transform_fingerprint=carried,
        source_sha256=str(payload.get("source_sha256") or ""),
        deck_digest=str(payload.get("deck_digest") or ""),
        evidence=dict((payload.get("finite_history") or {}).get("evidence")
                      or {}))
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
    # Compared in RECORDS, which is what the exporter's --increments bounds:
    # one increment of a C3D8 is eight probe records, one per integration
    # point. Comparing the increment count against the record count read a
    # complete 35-increment J2 history as 35 of 280 and refused it. The
    # exporter now writes both counts; the older key is the fallback for
    # fixtures frozen before it did.
    asked = frozen.get("records_requested", frozen.get("increments_requested"))
    got = frozen.get("records_carried", frozen.get("increments_carried"))
    if isinstance(asked, int) and isinstance(got, int) and got < asked:
        raise FixtureError(
            f"{path.name} carries {got} record(s) of the {asked} its "
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


def load_all(directory: Path, *,
             fingerprint: Optional[str] = CURRENT_TRANSFORM_FINGERPRINT
             ) -> list:
    """Every fixture in a directory, in a stable order."""
    return [load(path, fingerprint=fingerprint)
            for path in sorted(Path(directory).glob("*.json"))]


#: The two readings of "what is this UMAT's DDSDDE the derivative of".
#:
#: ``material``    Dsigma = D Denexpressed directly: the plain small-strain
#:                 reading, and what a routine that computes sigma = C:eps
#:                 returns even when the deck says NLGEOM=YES.
#: ``jaumann``     Dsigma = D Deps - sigma tr(Deps): Abaqus's finite-strain
#:                 material Jacobian, which is the tangent of the Jaumann rate
#:                 of KIRCHHOFF stress divided by J, so the Cauchy increment
#:                 carries the extra convective term.
TANGENT_READINGS = ("material", "jaumann")


@dataclass(frozen=True)
class TangentConvention:
    """Which reading of DDSDDE this fixture's own numbers support.

    Measured rather than inferred from the deck's ``NLGEOM``, because it is a
    property of the ROUTINE and not of the step: measured on the frozen set,
    ``irfancn__Abaqus-UMAT-elastic`` runs under ``NLGEOM=YES`` and its tangent
    satisfies the plain reading to 3.4e-16 while the Jaumann reading misses by
    5.5e-03 -- it is a small-strain routine handed a logarithmic strain --
    whereas ``AlexanderJFDR`` NeoHookean under the same flag is the other way
    round, 8.7e-03 against 6.9e-08. An assembler that picks by the flag builds
    the wrong stiffness for one of them.
    """

    reading: Optional[str]
    material_error: float
    jaumann_error: float
    separation: float
    increments_used: int
    volumetric: float
    why: str
    #: How much larger the tangent's own prediction is than the stress change
    #: it is supposed to predict. Far above one means the strain the fixture
    #: carries is not the strain this tangent is the derivative with respect
    #: to -- a growth or transformation model splits the total strain with an
    #: internal state and the stress answers only the elastic part.
    prediction_overshoot: float = 1.0

    @property
    def decided(self) -> bool:
        return self.reading is not None


def tangent_convention(fixture: VerifiedFixture, *,
                       tolerance: float = 1e-5) -> TangentConvention:
    """Which reading of the fixture's tangent differentiates its stress.

    Both readings are evaluated at the MID-POINT of each increment -- the
    average of the two reported tangents against the chord of the two reported
    stresses -- because a tangent at the end of an increment and a chord
    across it differ at first order in the increment for any nonlinear
    material, and that difference would be mistaken for a wrong tangent.

    A window in which the strain increment is isochoric cannot separate the
    two readings at all: the term they differ by is ``sigma tr(Deps)`` and the
    trace is zero. That is reported as undecided, not as agreement.
    """
    errors = {name: 0.0 for name in TANGENT_READINGS}
    used, volumetric, overshoot = 0, 0.0, 1.0
    for previous, current in zip(fixture.converted, fixture.converted[1:]):
        if current.tangent is None or previous.tangent is None:
            continue
        if current.dstrain.size != fixture.ntens:
            continue
        change = np.asarray(current.stress) - np.asarray(previous.stress)
        scale = float(np.max(np.abs(change)))
        if scale <= 0.0:
            continue
        dstrain = np.asarray(current.dstrain, dtype=float)
        middle_tangent = 0.5 * (current.tangent + previous.tangent)
        middle_stress = 0.5 * (np.asarray(previous.stress)
                               + np.asarray(current.stress))
        trace = float(np.sum(dstrain[:fixture.ndi])) if fixture.ndi else 0.0
        volumetric = max(volumetric, abs(trace))
        predicted = {
            "material": middle_tangent @ dstrain,
            "jaumann": middle_tangent @ dstrain - middle_stress * trace,
        }
        for name, value in predicted.items():
            errors[name] = max(errors[name],
                               float(np.max(np.abs(value - change))) / scale)
        overshoot = max(overshoot,
                        float(np.max(np.abs(predicted["material"]))) / scale)
        used += 1

    if not used:
        return TangentConvention(
            None, float("nan"), float("nan"), float("nan"), 0, 0.0,
            "no two increments of this window carry a tangent and a moving "
            "stress, so neither reading can be evaluated")
    material, jaumann = errors["material"], errors["jaumann"]
    if volumetric == 0.0:
        return TangentConvention(
            None, material, jaumann, 1.0, used, volumetric,
            "the window's strain increments are isochoric, and the two "
            "readings differ by sigma tr(Deps): with a zero trace they are the "
            "same formula and this window cannot tell them apart",
            prediction_overshoot=overshoot)
    best, other = ("material", "jaumann") if material <= jaumann \
        else ("jaumann", "material")
    separation = (errors[other] / errors[best]) if errors[best] > 0 else float("inf")
    if errors[best] > tolerance and overshoot > 10.0:
        return TangentConvention(
            None, material, jaumann, separation, used, volumetric,
            f"the tangent's own prediction is {overshoot:.3g} times the "
            f"stress change it is supposed to predict, so the strain this "
            f"fixture carries is not the strain this tangent differentiates "
            f"with respect to: a growth or transformation model splits the "
            f"total strain with an internal state and the stress answers only "
            f"the elastic part, which the fixture does not carry. Nothing "
            f"here says the tangent is wrong",
            prediction_overshoot=overshoot)
    if errors[best] > tolerance:
        return TangentConvention(
            None, material, jaumann, separation, used, volumetric,
            f"neither reading differentiates the reported stress: the plain "
            f"reading misses by {material:.3e} and the Jaumann one by "
            f"{jaumann:.3e}, both above {tolerance:.0e}",
            prediction_overshoot=overshoot)
    return TangentConvention(
        best, material, jaumann, separation, used, volumetric,
        f"the {best} reading holds to {errors[best]:.3e} over {used} "
        f"increment(s) and the {other} one misses by {errors[other]:.3e}, a "
        f"separation of {separation:.3g}", prediction_overshoot=overshoot)


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
