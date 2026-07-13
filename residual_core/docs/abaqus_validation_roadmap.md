# Abaqus Validation Roadmap

What is verified **offline today**, what needs **Abaqus** next, and the concrete
order to close each gap. Nothing here is a support claim beyond what is marked
verified.

## Legend

✅ offline-verified · 🟡 offline-capable once an export exists · ⏳ needs Abaqus
(+ Intel ifort for user subroutines) · ⛔ not implemented / license-limited

## Stage 0 — offline, done (no Abaqus)

| Item | State |
|---|---|
| Stress-driven C3D8 assembly == verified kernel (bitwise) | ✅ |
| Elastic material replay == kernel + FD tangent | ✅ |
| Truss `EA/L`, beam `PL³/3EI`, mixed dispatch | ✅ |
| Requirements engine, inspector, neutral IO, public API | ✅ |
| Validation scripts skip cleanly without Abaqus | ✅ |

## Stage 1 — Abaqus stress-driven cross-check (Mode 1)

Goal: confirm the framework's assembled reactions match Abaqus `RF` on a real,
**spatially-varying** stress state (validates integration-point ordering).

1. ⏳ `run_abaqus_validation.py` — run the internal C3D8 CP job.
2. ⏳ `extract_odb_fields.py` — export `S`, `RF`, `U`.
3. 🟡 `compare_residuals.py` — assemble offline, compare free-DOF residual ~ 0 and
   reactions vs `RF` (sign-aware). *This step is offline once the export exists.*

Unblocks: every external UMAT case (viscoelastic, J2, damage) via the same path.

## Stage 2 — external UMAT stress-driven (Mode 1)

For each permissive UMAT candidate (viscoelasticity, J2 plasticity, Mazars* ):
4. author a single-element `.inp` (some repos ship none), ⏳ run + export, 🟡 assemble.
Status flips `planned → implemented/offline-verified` per case once its export
exists. (* Mazars: resolve license caveat first.)

## Stage 3 — UMAT material replay (Mode 2)

5. ⏳ per-increment export (history), then compiled-UMAT replay (needs ifort) or an
   **independent** Python material replay (`compare_umat_replay.py`).
6. History-dependent laws: replay increment by increment (enforced).

## Stage 4 — UEL direct residual (Mode 3)

7. ⏳ instrument a permissive UEL (Elasticity/Hyperelasticity) to dump
   COORDS/U/PROPS/SVARS/RHS/AMATRX (`docs/uel_validation_plan.md`).
8. 🟡 `compare_uel_rhs.py` — sign convention, DOF order, AMATRX; then assemble.
   Flips `adapter-skeleton → verified` for `uel_direct`.

## Stage 5 — new backends (⛔ today)

| Backend | Blocker | Reference |
|---|---|---|
| Cohesive (Mode 2) | no separation-kinematics backend | GPL CZM (reference only) |
| Coupled-field (Mode 5) | no coupled DOF weak form | hydrogel UEL (custom license) |
| Shell (Mode 4) | contract only; no permissive Abaqus UEL | LS-Dyna shell (reference) |
| 2D / Cosserat UEL | 2D + micro-rotation DOF wiring | jgomezc1 (MIT) |

## Environment note

Abaqus is **not installed** in the current environment. Every ⏳ item is delivered
**ready-to-run and clearly marked pending**, never claimed as passing. The offline
suite (Stage 0) is fully green independent of Abaqus.
