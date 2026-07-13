"""Requirements engine — data minimization.

Answers three questions for a requested assembly *mode*, over a physics-free
"availability" snapshot (what the user has actually provided):

    What do we have?
    What do we need?
    What is the minimum next thing the user must provide?

It reports the **minimum missing item**, not a giant generic checklist. The modes
and their minimum inputs mirror ``docs/minimal_input_contract.md``. Nothing here
knows any mechanics — it reasons over declared input keys only.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class InputItem:
    key: str
    label: str
    recommendation: str = ""
    short_label: str = ""      # concise, user-facing status line (audit item 2)

    def display(self) -> str:
        return self.short_label or self.label


# Minimum inputs per mode (see docs/minimal_input_contract.md). The element field
# that enters the weak form is mode-dependent and intentionally generic: for
# solids it is integration-point stress, for beams/shells it is section
# resultants, for thermal it is heat flux — the engine does not care which.
MODE_REQUIREMENTS: Dict[str, List[InputItem]] = {
    "stress-driven": [
        InputItem("mesh", "mesh / connectivity / coordinates",
                  short_label="mesh"),
        InputItem("dof_field", "DOF solution field (e.g. displacement U, or "
                  "displacement+rotation for beams/shells, temperature for thermal)",
                  short_label="solution field (U / U+rotation / T)"),
        InputItem("element_field", "element field entering the weak form",
                  "integration-point stress field S, or an ODB/CSV export "
                  "(section resultants N/M/Q for beams/shells; heat flux for thermal)",
                  short_label="stress / resultant field"),
    ],
    "material-replay": [
        InputItem("mesh", "mesh / connectivity / coordinates",
                  short_label="mesh"),
        InputItem("solution_history", "solution history (not just the final step)",
                  "the solution history (the increment sequence, not just the "
                  "final step; history-dependent materials need it)",
                  short_label="solution history"),
        InputItem("material_model", "material model or adapter",
                  "a UMAT source, a built-in material, or a Python law",
                  short_label="material model / adapter"),
        InputItem("material_parameters", "material parameters (PROPS)",
                  short_label="material parameters (PROPS)"),
        InputItem("state_prev", "previous state variables (STATEV) or an initial state",
                  short_label="previous state / STATEV"),
        InputItem("time_increments", "time increments (dtime per step)",
                  short_label="time increments"),
    ],
    "direct-residual": [
        InputItem("mesh", "mesh / connectivity / coordinates",
                  short_label="mesh"),
        InputItem("element_dofs", "element DOF layout",
                  short_label="element DOF layout"),
        InputItem("uel_routine", "callable residual routine or UEL-like adapter",
                  "a routine returning (RHS, AMATRX, SVARS_new)",
                  short_label="callable UEL adapter"),
    ],
    "formulation": [
        InputItem("formulation_backend", "a formulation backend implementing the "
                  "base interface", "a registered Formulation subclass "
                  "(implement and register one)",
                  short_label="formulation backend"),
        InputItem("section_properties", "section / constitutive properties for the "
                  "formulation (e.g. E, A, I for a beam)",
                  "section properties (e.g. E, A for a bar; E, A, I for a beam)",
                  short_label="section properties (E, A, I)"),
    ],
}

# friendly aliases
_MODE_ALIASES = {
    "stress_driven": "stress-driven", "stressdriven": "stress-driven",
    "material_replay": "material-replay", "replay": "material-replay",
    "uel": "direct-residual", "uel-direct": "direct-residual",
    "direct": "direct-residual", "manufactured": "formulation",
    "manual": "formulation",
}


def canonical_mode(mode: str) -> str:
    m = (mode or "").strip().lower()
    return _MODE_ALIASES.get(m, m)


@dataclass
class RequirementsReport:
    mode: str
    have: List[InputItem] = field(default_factory=list)
    missing: List[InputItem] = field(default_factory=list)

    @property
    def runnable(self) -> bool:
        return not self.missing

    @property
    def minimum_next(self) -> Optional[InputItem]:
        return self.missing[0] if self.missing else None

    def render(self) -> str:
        """User-facing report. Shows availability and, when blocked, the single
        minimum missing input — never a giant generic checklist."""
        missing_keys = {it.key for it in self.missing}
        ordered = MODE_REQUIREMENTS.get(self.mode, self.have + self.missing)
        if self.runnable:
            L = ["Ready to assemble in %s mode." % self.mode, "Available:"]
            for it in ordered:
                L.append("  %s: yes" % it.display())
            L.append("")
            L.append("All minimum inputs are present.")
            return "\n".join(L)

        L = ["Cannot assemble in %s mode." % self.mode, "Available:"]
        for it in ordered:
            L.append("  %s: %s" % (it.display(), "no" if it.key in missing_keys else "yes"))
        nxt = self.minimum_next
        L.append("")
        L.append("Minimum missing input:")
        L.append("  provide %s." % (nxt.recommendation or nxt.label))
        return "\n".join(L)


def known_modes() -> List[str]:
    return sorted(MODE_REQUIREMENTS)


def evaluate_requirements(mode: str, available: Dict[str, bool]) -> RequirementsReport:
    """Given a mode and a {requirement_key -> bool} availability map, split the
    mode's minimum inputs into have / missing (missing keys default to absent)."""
    m = canonical_mode(mode)
    if m not in MODE_REQUIREMENTS:
        raise KeyError("unknown mode %r; known: %s" % (mode, ", ".join(known_modes())))
    report = RequirementsReport(mode=m)
    for item in MODE_REQUIREMENTS[m]:
        if bool(available.get(item.key, False)):
            report.have.append(item)
        else:
            report.missing.append(item)
    return report
