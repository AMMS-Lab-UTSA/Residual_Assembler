"""Which stage a bad residual came from, asked in the order the stages run.

A residual that is wrong has five places it can be wrong, and they need five
different people to fix them:

1. **the UMAT** -- the original routine returned something that is not a
   number, or a stress that does not move when the strain does;
2. **the transformation** -- the converted build computes a different stress
   from the original, so the two are not the same model and no derivative
   taken from one applies to the other;
3. **the constitutive derivative** -- the tangent the converted build reports
   is not the derivative of the stress it reports;
4. **the mapping** -- the tangent and the stress are right and the assembler
   is reading them in a different convention: a transposed Voigt order, or
   tensorial shear where the UMAT meant engineering;
5. **the residual assembly** -- everything handed in is right and the
   integration of ``B^T sigma`` is not.

"The residual is wrong by 12%" tells nobody which of those to open. Each check
below answers exactly one of them, in the order they run, and stops at the
first that fails -- because a stress that is NaN makes every later check
meaningless and reporting all five would bury the one that matters.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np

#: The stages, in the order they run. A failure is attributed to the first one
#: that does not hold, because everything after it is downstream of it.
STAGES = ("umat", "transformation", "constitutive_derivative", "mapping",
          "residual_assembly")

#: Who fixes each stage, in one line. Carried in the failure message because
#: the point of attributing a failure is to route it.
OWNER = {
    "umat": "the author's own routine, as published",
    "transformation": "the OTI source transformation",
    "constitutive_derivative": "the derivative the converted build extracts",
    "mapping": "the convention this assembler reads the UMAT's arrays in",
    "residual_assembly": "the integration of B^T sigma in this repository",
}


@dataclass
class Finding:
    """One stage's verdict, and the numbers behind it."""

    stage: str
    ok: bool
    detail: str = ""
    measured: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {"stage": self.stage, "ok": self.ok, "detail": self.detail,
                "measured": dict(self.measured)}

    def message(self) -> str:
        if self.ok:
            return f"{self.stage}: holds. {self.detail}"
        return (f"{self.stage} is where this goes wrong -- {OWNER[self.stage]}."
                f"\n  {self.detail}"
                + ("\n  measured: " + ", ".join(
                    f"{name}={value!r}" for name, value in self.measured.items())
                   if self.measured else ""))


@dataclass
class Diagnosis:
    """Every stage checked, and the first that failed."""

    findings: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return all(finding.ok for finding in self.findings)

    @property
    def blame(self) -> Optional[str]:
        for finding in self.findings:
            if not finding.ok:
                return finding.stage
        return None

    def report(self) -> str:
        lines = [finding.message() for finding in self.findings]
        if self.ok:
            return ("every stage holds: the UMAT computed numbers, the "
                    "conversion agrees with it, the tangent is the derivative "
                    "of the stress, the conventions match, and the assembly "
                    "integrates them.\n  " + "\n  ".join(lines))
        return (f"the failure is in {self.blame} -- {OWNER[self.blame]}.\n"
                + "\n".join(lines))

    def as_dict(self) -> dict:
        return {"ok": self.ok, "blame": self.blame,
                "findings": [f.as_dict() for f in self.findings]}


def _finite(values) -> bool:
    return bool(np.all(np.isfinite(np.asarray(values, dtype=float))))


def _relative(a, b) -> float:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    scale = max(float(np.max(np.abs(a))), float(np.max(np.abs(b))), 0.0)
    if scale == 0.0:
        return 0.0
    return float(np.max(np.abs(a - b)) / scale)


def diagnose(fixture, *, assemble: Callable, reference: Callable,
             primal_tolerance: float = 1e-9,
             derivative_tolerance: float = 1e-5,
             assembly_tolerance: float = 1e-10) -> Diagnosis:
    """Run every stage's check against one verified fixture.

    ``assemble(stress_at_ip)`` returns this repository's element internal
    force. ``reference(stress_at_ip)`` returns the same quantity computed some
    other way -- by hand, or by a second implementation -- so that stage five
    is checked against something rather than against itself.
    """
    findings: list = []

    # 1. Did the author's routine compute numbers, and do they respond?
    stresses = [record.stress for record in fixture.original]
    strains = [record.strain for record in fixture.original]
    non_finite = [record.increment for record in fixture.original
                  if not _finite(record.stress) or not _finite(record.state)]
    if non_finite:
        findings.append(Finding(
            "umat", False,
            f"the ORIGINAL routine returned values that are not numbers at "
            f"increment(s) {non_finite}. Nothing downstream of a NaN is a "
            f"measurement.",
            {"increments": non_finite}))
        return Diagnosis(findings)
    moved = _relative(stresses[0], stresses[-1]) if len(stresses) > 1 else 0.0
    if len(stresses) > 1 and moved == 0.0 and _relative(strains[0], strains[-1]) > 0:
        findings.append(Finding(
            "umat", False,
            "the strain changed over the history and the stress did not, so "
            "the routine is not responding to what it is given.",
            {"strain_change": _relative(strains[0], strains[-1])}))
        return Diagnosis(findings)
    findings.append(Finding(
        "umat", True,
        f"{len(stresses)} increments, all finite, stress responding to strain",
        {"stress_change": moved}))

    # 2. Does the converted build compute the same stress?
    worst = max((_relative(a.stress, b.stress)
                 for a, b in zip(fixture.original, fixture.converted)),
                default=0.0)
    worst_state = max((_relative(a.state, b.state)
                       for a, b in zip(fixture.original, fixture.converted)
                       if a.state.size), default=0.0)
    if worst > primal_tolerance or worst_state > primal_tolerance:
        findings.append(Finding(
            "transformation", False,
            f"the converted build's stress differs from the original's by "
            f"{worst:.3e} and its state by {worst_state:.3e}. They are not the "
            f"same model, so no derivative from one describes the other.",
            {"stress": worst, "state": worst_state}))
        return Diagnosis(findings)
    findings.append(Finding(
        "transformation", True,
        f"original and converted agree to {worst:.3e} in stress and "
        f"{worst_state:.3e} in state over {fixture.increments()} increments"))

    # 3. Is the reported tangent the derivative of the reported stress?
    worst_derivative, at = _derivative_check(fixture)
    if worst_derivative is None:
        findings.append(Finding(
            "constitutive_derivative", True,
            "the history has no two increments whose strain increment is "
            "large enough to difference; nothing is claimed either way"))
    elif worst_derivative > derivative_tolerance:
        findings.append(Finding(
            "constitutive_derivative", False,
            f"the tangent the converted build reports does not predict the "
            f"stress increment it reports: {worst_derivative:.3e} relative at "
            f"increment {at}. A tangent that does not differentiate its own "
            f"stress is not a tangent.",
            {"relative": worst_derivative, "increment": at}))
        return Diagnosis(findings)
    else:
        findings.append(Finding(
            "constitutive_derivative", True,
            f"D dstrain predicts the stress increment to "
            f"{worst_derivative:.3e} over the history"))

    # 4. Do the two sides read the arrays the same way?
    from residual_core.materials.verified_fixture import check_conventions

    problems = check_conventions(fixture)
    if problems:
        findings.append(Finding(
            "mapping", False,
            "the fixture's own numbers contradict the convention it states: "
            + "; ".join(problems)))
        return Diagnosis(findings)
    findings.append(Finding(
        "mapping", True,
        f"{fixture.ntens} components as {fixture.ndi} direct and "
        f"{fixture.nshr} engineering shear, and the tangent is shaped and "
        f"signed accordingly"))

    # 5. Does the assembly integrate them?
    stress_at_ip = np.tile(fixture.original[-1].stress, (8, 1))
    mine = np.asarray(assemble(stress_at_ip), dtype=float)
    theirs = np.asarray(reference(stress_at_ip), dtype=float)
    difference = _relative(mine, theirs)
    if difference > assembly_tolerance:
        findings.append(Finding(
            "residual_assembly", False,
            f"the assembled internal force differs from the independent "
            f"reference by {difference:.3e}. Everything handed to it agreed, "
            f"so this is the integration.",
            {"relative": difference}))
        return Diagnosis(findings)
    findings.append(Finding(
        "residual_assembly", True,
        f"B^T sigma integrates to the independent reference within "
        f"{difference:.3e}"))
    return Diagnosis(findings)


def _derivative_check(fixture) -> tuple:
    """Does ``D dstrain`` predict the stress increment the same build reported?

    A weaker question than the one the UMAT pipeline asks -- it perturbs the
    routine and this only reads what it wrote -- and a different one: it asks
    whether the tangent and the stress in THIS fixture belong together. That
    is the question the assembler needs answered, because it uses both.
    """
    worst, at = None, None
    for previous, current in zip(fixture.converted, fixture.converted[1:]):
        if current.tangent is None or current.dstrain.size != fixture.ntens:
            continue
        change = np.asarray(current.stress) - np.asarray(previous.stress)
        scale = float(np.max(np.abs(change)))
        if scale <= 0.0:
            continue
        predicted = current.tangent @ np.asarray(current.dstrain, dtype=float)
        relative = _relative(predicted, change)
        if worst is None or relative > worst:
            worst, at = relative, current.increment
    return worst, at
