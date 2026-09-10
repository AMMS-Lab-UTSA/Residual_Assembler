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

The other thing it carries is the CONVENTION. A tangent is not a number until
somebody says what its indices mean, and the two sides of this bridge have to
agree about Voigt ordering and about which off-diagonal entries are
engineering shear. The fixture states it; :func:`check_conventions` checks the
numbers are consistent with what it states, so a mapping error is caught where
it happens rather than as a residual that is wrong by a factor of two.
"""
from __future__ import annotations

import json
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
    for index, record in enumerate(fixture.original):
        if record.stress.size != ntens:
            raise FixtureError(
                f"{path.name} increment {record.increment}: the stress has "
                f"{record.stress.size} components and the fixture states "
                f"NTENS={ntens}. One of the two is wrong and neither may be "
                f"assumed.")
    return fixture


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
