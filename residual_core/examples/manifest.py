"""Formulation-aware classifier / regenerator for the example manifest.

The Step-1 manifest (``residual_core/example_manifest.json``) inventoried every
Abaqus ``.inp`` under ``sources/`` and grouped them with the original
crystal-plasticity taxonomy (old ``group`` A/B/C/D). That taxonomy was too
narrow for the refactored, formulation-agnostic core: it keyed on *which CP
milestone* a file belonged to rather than *which residual backend* would
assemble it.

This module re-classifies every recorded example **by formulation / residual
backend** and regenerates ``example_manifest.json`` + ``example_manifest.md`` in
place, without re-parsing the raw ``.inp`` files (it reuses the element types,
user-material / user-element flags, constant / state counts, and include
resolution already recorded in Step-1). No file under ``sources/`` is touched.

Primary grouping — ``formulation_group`` (exactly one per example, in this
precedence):

  F. Reference-only due to license  -- ANY file under ``sources/copyleft/`` or
     ``sources/license-unknown/`` (hard licensing gate, regardless of technical
     shape). Its technical shape is still recorded in ``technical_group``.
  C. UEL / direct-residual          -- uses ``*User Element`` or cohesive
     elements (``U1``, ``COH*``).
  D. Coupled-field                  -- coupled temperature-displacement /
     thermal-mechanical / pore-pressure / piezo (element type ending in ``T``,
     ``*Coupled Temperature-displacement``, ...).
  A. Standard solid element + UMAT  -- standard continuum solid element with a
     ``*User Material``; not C/D.
  B. Standard solid element + built-in material -- standard continuum solid,
     built-in material (no ``*User Material``); not C/D.
  E. Unsupported element type       -- a standard analysis whose element type is
     not a continuum solid we can place into A/B, or a mesh/BC/material
     *fragment* / meshless / special element. Used only when A/B don't fit.

``technical_group`` records the same A-E *shape* IGNORING the licensing gate, so
reference-only (F) files still show what they technically are. The pre-existing
Step-1 ``technical_group`` value is preserved as ``technical_group_legacy``.

Public API
----------
    classify(record) -> dict          # the new formulation-aware fields
    enrich(record)   -> dict          # record + new fields (non-destructive)
    generate(existing_json_path, out_json, out_md) -> dict   # regenerate both

Run ``python residual_core/examples/manifest.py`` to regenerate
``residual_core/example_manifest.json`` and ``residual_core/example_manifest.md``
in place (and print the F-invariant + priority-file checks).
"""

from __future__ import annotations

import json
import os
import re
from collections import OrderedDict

# --------------------------------------------------------------------------- #
# Locations (this file lives in residual_core/examples/)
# --------------------------------------------------------------------------- #
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_CORE_DIR = os.path.dirname(_THIS_DIR)                       # residual_core/
DEFAULT_JSON = os.path.join(_CORE_DIR, "example_manifest.json")
DEFAULT_MD = os.path.join(_CORE_DIR, "example_manifest.md")

# --------------------------------------------------------------------------- #
# Vocabulary
# --------------------------------------------------------------------------- #
LICENSE_REFERENCE_ONLY = ("copyleft", "license-unknown")

FORMULATION_GROUPS = ("A", "B", "C", "D", "E", "F")
FORMULATION_GROUP_LABELS = OrderedDict([
    ("F", "Reference-only due to license"),
    ("C", "UEL / direct-residual"),
    ("D", "Coupled-field"),
    ("A", "Standard solid elements + UMAT"),
    ("B", "Standard solid elements + built-in material"),
    ("E", "Unsupported element type / fragment"),
])

# Verification levels (fallback ladder from the task; docs/verification_strategy.md
# is not present yet). 0 parser, 1 zero-field, 2 rigid-body, 3 patch,
# 4 FD-tangent, 5 solver-vs-Abaqus, 6 material-replay, 7 full-replay.
VERIFICATION_LEVEL_NAMES = OrderedDict([
    (0, "parser"),
    (1, "zero-field"),
    (2, "rigid-body"),
    (3, "patch"),
    (4, "FD-tangent"),
    (5, "solver-vs-Abaqus"),
    (6, "material-replay"),
    (7, "full-replay"),
])

# The 4 priority C3D8 crystal-plasticity targets (must land in group A /
# solid_c3d8_finite_strain / support_status "supported").
PRIORITY_TARGETS = (
    "sources/permissive/ngrilli_Oxford_Crystal_Plasticity/ExampleInputFiles/HCPnoTwin/Compression111.inp",
    "sources/permissive/ngrilli_Oxford_Crystal_Plasticity/ExampleInputFiles/DiscreteTwin/Job-1.inp",
    "sources/permissive/ngrilli_Oxford_Crystal_Plasticity/ExampleInputFiles/HomogeneousTwin/Job-1.inp",
    "sources/permissive/ngrilli_Oxford_Crystal_Plasticity/InterfaceElementGeneration/PyCiGen/Job-1.inp",
)


# --------------------------------------------------------------------------- #
# Element-type predicates
# --------------------------------------------------------------------------- #
def is_user_element_type(t):
    """User element (Abaqus UEL) or cohesive family: U<n>, COH*."""
    t = (t or "").upper()
    return bool(re.match(r"^U\d+$", t)) or t.startswith("COH")


def is_coupled_type(t):
    """Continuum element carrying a temperature DOF (thermo-mechanical), e.g.
    C3D8T, CPE4T, CAX8RT -> name ends in 'T' on a continuum family."""
    t = (t or "").upper()
    return t.endswith("T") and bool(re.match(r"^(C3D|CPE|CPS|CAX|CPEG)", t))


def is_continuum_solid(t):
    """Standard continuum solid element (3D or 2D), excluding coupled ('..T')."""
    t = (t or "").upper()
    return bool(re.match(r"^(C3D|CPE|CPS|CAX|CPEG)\d", t)) and not t.endswith("T")


def is_3d_solid(t):
    return (t or "").upper().startswith("C3D")


def is_2d_solid(t):
    t = (t or "").upper()
    return bool(re.match(r"^(CPE|CPS|CAX|CPEG)", t))


# --------------------------------------------------------------------------- #
# Effective (include-resolved) view of a record
# --------------------------------------------------------------------------- #
def effective_fields(record):
    """Return the effective element types / user-material / user-element flags,
    folding in ``resolved_from_includes`` for multi-file masters. Reuses only
    fields already present in the Step-1 record (no raw re-parsing)."""
    et = list(record.get("element_types") or [])
    uses_um = bool(record.get("uses_user_material"))
    uses_ue = bool(record.get("uses_user_element"))

    resolved = record.get("resolved_from_includes") or {}
    if resolved:
        if not et and resolved.get("element_types"):
            et = list(resolved.get("element_types") or [])
        uses_um = uses_um or bool(resolved.get("uses_user_material"))
        uses_ue = uses_ue or bool(resolved.get("uses_user_element"))

    # user_element_types / cohesive_element_types are additional evidence
    if record.get("user_element_types"):
        uses_ue = True
    return et, uses_um, uses_ue


def _has_builtin_material(record):
    if (record.get("n_materials") or 0) > 0:
        return True
    if record.get("materials"):
        return True
    resolved = record.get("resolved_from_includes") or {}
    if resolved.get("n_material_constants"):
        return True
    return False


def _is_coupled(record, etypes):
    if any(is_coupled_type(t) for t in etypes):
        return True
    notes = (record.get("notes") or "").lower()
    return ("coupled temp" in notes or "pore pressure" in notes or
            "piezo" in notes or "thermal-mechanical" in notes or
            "thermomechanical" in notes)


# --------------------------------------------------------------------------- #
# Core classification
# --------------------------------------------------------------------------- #
def technical_group(record):
    """A-E technical shape IGNORING the licensing gate.

    Precedence: C (UEL) > D (coupled) > A/B (standard continuum solid,
    UMAT vs built-in) > E (fragment / unsupported / special)."""
    etypes, uses_um, uses_ue = effective_fields(record)
    cohesive = record.get("cohesive_element_types") or []
    user_elems = record.get("user_element_types") or []

    if uses_ue or user_elems or cohesive or any(is_user_element_type(t) for t in etypes):
        return "C"
    if _is_coupled(record, etypes):
        return "D"
    if any(is_continuum_solid(t) for t in etypes):
        return "A" if uses_um else "B"
    return "E"


def formulation_group(record):
    """Primary grouping. Licensing gate (F) dominates the technical shape."""
    if record.get("license_tier") in LICENSE_REFERENCE_ONLY:
        return "F"
    return technical_group(record)


def _c3d8_only(etypes):
    s = {(t or "").upper() for t in etypes}
    return s == {"C3D8"}


def dofs_for(record, tgroup, etypes):
    """Per-node DOF list. Solids -> UX/UY/UZ (2D -> UX/UY); coupled adds T;
    UEL and fragments -> 'unknown' (unless trivially inferable)."""
    if tgroup == "C":
        return "unknown"
    if tgroup == "E":
        return "unknown"
    base = None
    if any(is_3d_solid(t) for t in etypes):
        base = ["UX", "UY", "UZ"]
    elif any(is_2d_solid(t) for t in etypes):
        base = ["UX", "UY"]
    if base is None:
        return "unknown"
    if tgroup == "D":
        return base + ["T"]
    return base


def material_type_for(record, tgroup, uses_um):
    if uses_um:
        return "umat"
    if tgroup == "C":
        return "uel-internal"
    if _has_builtin_material(record):
        return "builtin"
    return "none"


def residual_backend_for(record, tgroup, etypes, nlgeom):
    """Which implemented framework backend would assemble this example
    (technical, license-agnostic)."""
    if tgroup == "C":
        return "uel_direct"
    if tgroup == "D":
        return "coupled (not implemented)"
    if tgroup == "E":
        return "none (unsupported element)"
    # A / B : standard continuum solid
    if _c3d8_only(etypes):
        if tgroup == "A":
            return ("solid_c3d8_finite_strain" if nlgeom
                    else "solid_c3d8_small_strain")
        # built-in-material C3D8: residual assembled from exported stress
        return "stress_driven_c3d8"
    # standard solid but not the C3D8 kernel (C3D8R/C3D20R/C3D4/CAX4/CPE*/...)
    return "none (unsupported element)"


def support_status_for(fgroup, tgroup, etypes):
    if fgroup == "F":
        return "reference-only"
    if tgroup in ("A", "B"):
        return "supported" if _c3d8_only(etypes) else "planned"
    if tgroup == "C":
        return "adapter-scaffolded"
    if tgroup == "D":
        return "not-implemented"
    return "not-implemented"          # E: fragment / special element


def required_verification_levels(fgroup, tgroup):
    """Subset of [0..7] that applies. Reference-only (F) collapses to parser
    only; otherwise keyed on the technical shape."""
    if fgroup == "F":
        return [0]
    if tgroup == "A":                     # UMAT solid
        return [0, 1, 2, 3, 4, 5, 6, 7]
    if tgroup == "B":                     # built-in-material solid
        return [0, 1, 2, 3, 4, 5]
    if tgroup == "C":                     # UEL
        return [0, 1, 4, 5, 7]
    if tgroup == "D":                     # coupled (not implemented)
        return [0]
    return [0]                            # E: fragment / special


def classify(record):
    """Return the new formulation-aware fields for one manifest record.

    Does not mutate ``record``. Reuses the Step-1 parsed fields (element types,
    user-material / user-element flags, counts, include resolution)."""
    etypes, uses_um, uses_ue = effective_fields(record)
    nlgeom = bool(record.get("nlgeom"))

    tgroup = technical_group(record)
    fgroup = formulation_group(record)

    return {
        "formulation_group": fgroup,
        "technical_group": tgroup,
        "dofs": dofs_for(record, tgroup, etypes),
        "material_type": material_type_for(record, tgroup, uses_um),
        "residual_backend": residual_backend_for(record, tgroup, etypes, nlgeom),
        "support_status": support_status_for(fgroup, tgroup, etypes),
        "required_verification_levels":
            required_verification_levels(fgroup, tgroup),
    }


def enrich(record):
    """Return a copy of ``record`` with the new fields merged in, preserving all
    existing data. The pre-existing ``technical_group`` (Step-1 legacy taxonomy)
    is retained as ``technical_group_legacy``. Idempotent across re-runs."""
    out = dict(record)
    if "technical_group_legacy" not in out and "technical_group" in out:
        out["technical_group_legacy"] = out["technical_group"]
    out.update(classify(record))
    return out


# --------------------------------------------------------------------------- #
# Regeneration
# --------------------------------------------------------------------------- #
def _counts(examples, key):
    c = {}
    for e in examples:
        c[e[key]] = c.get(e[key], 0) + 1
    return OrderedDict(sorted(c.items()))


def _ordered_groups(counts):
    """Formulation groups in the documentation order F, A, B, C, D, E."""
    order = ["F", "A", "B", "C", "D", "E"]
    return [g for g in order if counts.get(g)]


def generate(existing_json_path=DEFAULT_JSON, out_json=DEFAULT_JSON,
             out_md=DEFAULT_MD):
    """Read the existing manifest JSON, enrich every record with the
    formulation-aware classification, and (over)write both ``out_json`` and
    ``out_md``. Returns the regenerated top-level dict."""
    with open(existing_json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    examples = [enrich(e) for e in data.get("examples", [])]

    counts_form = _counts(examples, "formulation_group")
    counts_tech = _counts(examples, "technical_group")

    out = OrderedDict()
    # Preserve/refresh top-level metadata, keeping legacy summaries.
    out["generated_note"] = (
        "Formulation-aware manifest of every Abaqus .inp under "
        "sources/{permissive,copyleft,license-unknown}. Regenerated by "
        "residual_core/examples/manifest.py, which re-classifies the Step-1 "
        "records BY FORMULATION / RESIDUAL BACKEND (field formulation_group, "
        "A-F). The licensing gate F dominates: every file under copyleft/ or "
        "license-unknown/ is group F (reference-only) regardless of technical "
        "shape, and its shape is recorded separately in technical_group. All "
        "Step-1 fields are retained; the old Step-1 technical_group value is "
        "kept as technical_group_legacy and the old primary taxonomy as group. "
        "Nothing under sources/ was modified.")
    out["generated_date"] = data.get("generated_date", "2026-07-10")
    out["taxonomy"] = OrderedDict(
        (g, FORMULATION_GROUP_LABELS[g]) for g in ["F", "C", "D", "A", "B", "E"])
    out["verification_levels"] = OrderedDict(
        (str(k), v) for k, v in VERIFICATION_LEVEL_NAMES.items())
    out["total_inp_files"] = len(examples)
    out["counts_by_formulation_group"] = counts_form
    out["counts_by_technical_group"] = counts_tech
    # legacy summaries (kept, not lost) — recomputed from the per-example
    # legacy fields so this stays correct across repeated regenerations.
    if any("group" in e for e in examples):
        out["counts_by_group_legacy"] = _counts(
            [e for e in examples if "group" in e], "group")
    if any("technical_group_legacy" in e for e in examples):
        out["counts_by_technical_group_legacy"] = _counts(
            [e for e in examples if "technical_group_legacy" in e],
            "technical_group_legacy")
    out["priority_targets"] = list(PRIORITY_TARGETS)
    out["examples"] = examples

    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")

    _write_md(out, out_md)
    return out


# --------------------------------------------------------------------------- #
# Markdown rendering
# --------------------------------------------------------------------------- #
def _rel(path):
    return path


def _md_table_row(cols):
    return "| " + " | ".join(str(c) for c in cols) + " |"


def _write_md(out, out_md):
    examples = out["examples"]
    counts_form = out["counts_by_formulation_group"]
    total = out["total_inp_files"]

    by_group = {g: [] for g in FORMULATION_GROUPS}
    for e in examples:
        by_group[e["formulation_group"]].append(e)

    L = []
    L.append("# Example Manifest — classification by formulation / residual backend")
    L.append("")
    L.append(out["generated_note"])
    L.append("")
    L.append("Generated: %s. Total `.inp` files: **%d**." %
             (out["generated_date"], total))
    L.append("")

    # ---- summary counts table by formulation_group ----
    L.append("## Summary — counts by `formulation_group`")
    L.append("")
    L.append(_md_table_row(["Group", "Meaning", "Count", "Residual backend",
                            "Support status"]))
    L.append(_md_table_row(["---", "---", "---:", "---", "---"]))
    _backend_hint = {
        "F": "n/a (reference-only)",
        "C": "uel_direct",
        "D": "coupled (not implemented)",
        "A": "solid_c3d8_finite/small_strain",
        "B": "stress_driven_c3d8 / (planned)",
        "E": "none (unsupported element)",
    }
    _support_hint = {
        "F": "reference-only",
        "C": "adapter-scaffolded",
        "D": "not-implemented",
        "A": "supported (C3D8) / planned",
        "B": "supported (C3D8) / planned",
        "E": "not-implemented",
    }
    for g in _ordered_groups(counts_form):
        L.append(_md_table_row([g, FORMULATION_GROUP_LABELS[g],
                                counts_form.get(g, 0), _backend_hint[g],
                                _support_hint[g]]))
    L.append(_md_table_row(["**Total**", "", "**%d**" % total, "", ""]))
    L.append("")

    # ---- priority callout ----
    prio = [e for e in examples if e["path"] in set(PRIORITY_TARGETS)]
    L.append("## Priority targets — the 4 C3D8 crystal-plasticity files")
    L.append("")
    L.append("These are the first-milestone targets: standard continuum **C3D8** "
             "elements + crystal-plasticity **UMAT**, `nlgeom=YES` (finite "
             "strain). They MUST classify as group **A**, backend "
             "**`solid_c3d8_finite_strain`**, support **`supported`**.")
    L.append("")
    L.append(_md_table_row(["File", "formulation_group", "residual_backend",
                            "support_status"]))
    L.append(_md_table_row(["---", "---", "---", "---"]))
    order = {p: i for i, p in enumerate(PRIORITY_TARGETS)}
    for e in sorted(prio, key=lambda r: order.get(r["path"], 99)):
        L.append(_md_table_row(["`%s`" % _rel(e["path"]),
                                e["formulation_group"], e["residual_backend"],
                                e["support_status"]]))
    L.append("")
    L.append("> Note: **license-F files also carry a `technical_group`** — the "
             "same A-E shape computed while ignoring the licensing gate — so a "
             "reference-only file still shows what it technically is (e.g. a "
             "copyleft C3D8+UMAT job is `formulation_group=F`, "
             "`technical_group=A`). The pre-Step-1 taxonomy is preserved per "
             "record as `technical_group_legacy` and `group`.")
    L.append("")

    # ---- per-group sections ----
    L.append("## Per-group detail")
    L.append("")
    for g in _ordered_groups(counts_form):
        rows = by_group[g]
        L.append("### Group %s — %s (%d)" %
                 (g, FORMULATION_GROUP_LABELS[g], len(rows)))
        L.append("")
        if g == "F":
            L.append("Reference-only: under `sources/copyleft/` or "
                     "`sources/license-unknown/`. Never integrated or linked; "
                     "verification collapses to level 0 (parser). The "
                     "`technical_group` column shows the shape they would have "
                     "if licensing allowed.")
        L.append("")
        L.append(_md_table_row(["path", "technical_group", "residual_backend",
                                "material_type", "dofs", "support_status",
                                "verif. levels"]))
        L.append(_md_table_row(["---"] * 7))
        for e in sorted(rows, key=lambda r: r["path"]):
            dofs = e["dofs"]
            dofs_s = "/".join(dofs) if isinstance(dofs, list) else dofs
            lv = ",".join(str(i) for i in e["required_verification_levels"])
            L.append(_md_table_row(["`%s`" % _rel(e["path"]),
                                    e["technical_group"], e["residual_backend"],
                                    e["material_type"], dofs_s,
                                    e["support_status"], lv]))
        L.append("")

    # ---- verification level legend ----
    L.append("## Verification levels")
    L.append("")
    for k, v in VERIFICATION_LEVEL_NAMES.items():
        L.append("- **%d** — %s" % (k, v))
    L.append("")
    L.append("UMAT solids (A) require 0-7; built-in-material solids (B) 0-5; "
             "UEL (C) 0,1,4,5,7; coupled (D) and fragments (E) parser only; "
             "reference-only (F) parser only.")
    L.append("")

    with open(out_md, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))


# --------------------------------------------------------------------------- #
# Self-check helpers
# --------------------------------------------------------------------------- #
def check_invariants(out):
    """Verify: (1) every copyleft/license-unknown file is group F; (2) the 4
    priority files are A / solid_c3d8_finite_strain / supported. Returns a dict
    of results and prints PASS/FAIL lines."""
    examples = out["examples"]
    results = {}

    violations = [e["path"] for e in examples
                  if e.get("license_tier") in LICENSE_REFERENCE_ONLY
                  and e["formulation_group"] != "F"]
    results["f_invariant_ok"] = not violations
    results["f_invariant_violations"] = violations
    print("F-invariant (every copyleft/license-unknown file is group F): %s"
          % ("PASS" if not violations else "FAIL (%d)" % len(violations)))

    prio = {e["path"]: e for e in examples}
    ok_prio = True
    for p in PRIORITY_TARGETS:
        e = prio.get(p)
        good = (e is not None and e["formulation_group"] == "A"
                and e["residual_backend"] == "solid_c3d8_finite_strain"
                and e["support_status"] == "supported")
        ok_prio = ok_prio and good
    results["priority_ok"] = ok_prio
    print("Priority files A / solid_c3d8_finite_strain / supported: %s"
          % ("PASS" if ok_prio else "FAIL"))
    return results


def main():
    out = generate(DEFAULT_JSON, DEFAULT_JSON, DEFAULT_MD)
    print("Regenerated:")
    print("  %s" % DEFAULT_JSON)
    print("  %s" % DEFAULT_MD)
    print("total_inp_files: %d" % out["total_inp_files"])
    print("counts_by_formulation_group: %s"
          % json.dumps(out["counts_by_formulation_group"]))
    print("counts_by_technical_group:   %s"
          % json.dumps(out["counts_by_technical_group"]))
    check_invariants(out)


if __name__ == "__main__":
    main()
