# Verification Zoo

A **map of external material models and user elements** used to check the
framework's *formulation-agnostic* claim across many physics — not only crystal
plasticity. It is an honest landscape, not a support claim: most cases are
`planned`, `reference-only`, or `unsupported`, and say exactly what is missing.

Provenance and licenses live in [`sources/external_manifest.md`](../../sources/external_manifest.md),
[`sources/external_manifest.json`](../../sources/external_manifest.json), and
[`sources/LICENSE_AUDIT.md`](../../sources/LICENSE_AUDIT.md). **No external code is
vendored**; copyleft/unknown/different-solver sources are reference-only.

## Layout

```
tests/verification_zoo/
  README.md              (this file)
  manifest_summary.md    one-line status per case
  umat_cases/            standard-element + UMAT candidates (Class B)
  uel_cases/             UEL / direct-residual candidates (Class D)
  stress_driven_cases/   Mode-1 exported-field verification candidates (Class F)
  unsupported_cases/     no backend yet / license-limited (Class E, G, H)
    <case>/CASE.md        one verification card per case
```

## Residual-path classification (per card)

| Class | Meaning |
|---|---|
| A | standard element + built-in material |
| B | standard element + UMAT |
| C | standard element + VUMAT |
| D | UEL / direct residual |
| E | coupled-field formulation |
| F | stress-driven verification candidate |
| G | unsupported but useful later |
| H | reference-only due to license |

## Minimum verification path (per card)

| Mode | How |
|---|---|
| 1 stress-driven | assemble residual from exported solver fields (stress/flux/resultants/traction) |
| 2 material-replay | replay the UMAT-like update outside Abaqus, then assemble |
| 3 direct UEL | call/compare a UEL-like RHS/AMATRX |
| 4 native formulation | implement/verify the weak form directly (truss, beam, shell, thermal) |

## Status vocabulary

`implemented/offline-verified` · `implemented/Abaqus-pending` ·
`adapter-skeleton` · `planned` · `unsupported` · `reference-only/license-limited`

## How to use

- Offline, today: run the framework's own suite (see repository
  [`REVIEW_GUIDE.md`](../../REVIEW_GUIDE.md)). The zoo cards do **not** require
  Abaqus and do **not** run external code.
- With Abaqus later: follow [`docs/umat_validation_plan.md`](../../residual_core/docs/umat_validation_plan.md)
  and [`docs/uel_validation_plan.md`](../../residual_core/docs/uel_validation_plan.md),
  using the ready-to-run scripts under [`scripts/`](../../scripts) (they skip
  cleanly when Abaqus is absent).
