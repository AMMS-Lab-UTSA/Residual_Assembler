# License Audit — External Verification Candidates

Audit date: **2026-07-13**. Each external candidate's license was read from its
public GitHub repository (About panel / LICENSE file) on the review date. This
document records the determination and the **integration decision** for each.

## Decision rules

| License class | Examples | Integration decision |
|---|---|---|
| Permissive | MIT, BSD-2/3-Clause, Apache-2.0 | May be vendored **with attribution + license retention**. Not required for offline tests. |
| Copyleft | GPL-2.0/3.0, AGPL-3.0, LGPL | **Reference only.** Never copy/adapt/link/import into the framework. Re-derive from published equations if ever needed. |
| Custom / non-SPDX | bespoke LICENSE text | **Reference only** until the exact terms are reviewed and confirmed to permit reuse. |
| None / unclear | no LICENSE file | **All rights reserved** by default → reference only. |
| Different solver | LS-Dyna, ANSYS, COMSOL-only | Conceptual reference; not an Abaqus source. |

## Per-candidate determination

| Candidate | License (verified) | Class | Placement | Integrate? | Rationale |
|---|---|---|---|---|---|
| thealanjason/umat_finite_viscoelasticity | **MIT** | permissive | external_candidates | eligible (later) | Clear MIT; permissive; not needed for offline tests yet |
| sergiolucarini/UMAT4COMSOL | **MIT** | permissive | external_candidates | eligible (later) | MIT stated in README + About; J2 plasticity source |
| marioruiarruda/Mazars_UMAT | **BSD-3-Clause** ⚠ | permissive-with-caveat | external_candidates | **hold** | BSD-3 LICENSE, but extra `1-Permissions_&_Instructions.txt` — review author terms first |
| Sina-Taghizadeh/UMAT_Hyperelastic | **GPL-3.0** | copyleft | reference_only | **no** | GPL copyleft would virally relicense the framework |
| harshaa765/Bilinear-CZM-UMAT | **GPL-3.0** | copyleft | reference_only | **no** | GPL copyleft; cohesive backend not implemented anyway |
| bibekananda-datta/Abaqus-UEL-Elasticity | **BSD-3-Clause** | permissive | permissive/ (mirrored) | eligible | Already mirrored; clean BSD-3 |
| bibekananda-datta/Abaqus-UEL-Hyperelasticity | **BSD-3-Clause** | permissive | permissive/ (mirrored) | eligible | Already mirrored; clean BSD-3 |
| bibekananda-datta/Abaqus-UEL-Hydrogel | **custom (non-SPDX)** | unknown | reference_only | **no (yet)** | "View license" = bespoke terms; confirm before reuse |
| jgomezc1/ABAQUS-US | **MIT** | permissive | permissive/ (mirrored) | eligible | Already mirrored; MIT |
| CAEAssistant-Group/Abaqus-UEL-Subroutine | **MIT** (demo) | permissive-but-incomplete | reference_only | **no** | MIT, but intentionally partial teaser; not runnable |
| jfriedlein/usld·ushl_LS-Dyna_Fortran | unverified | different solver | reference_only | **no** | LS-Dyna, not Abaqus; conceptual reference only |

## Copyleft / restricted flags

- **GPL-3.0**: `UMAT_Hyperelastic`, `Bilinear-CZM-UMAT`. Isolated under
  `sources/reference_only/`. Must not be compiled, `include`d, linked, or copied
  into any framework code. If a hyperelastic or cohesive material is ever wanted,
  it must be **re-implemented independently** from the published equations.
- **Custom/non-SPDX**: `Abaqus-UEL-Hydrogel`. The LICENSE is a bespoke file, not a
  recognized SPDX identifier; treat as all-rights-reserved until the exact grant
  is confirmed.
- **BSD-3 with extra conditions**: `Mazars_UMAT` ships an author
  permissions/instructions file alongside the BSD-3 LICENSE. Resolve the apparent
  double-terms before any reuse.

## What is safe to build against

Only permissive, SPDX-clear sources with confirmed terms:
`Abaqus-UEL-Elasticity` (BSD-3), `Abaqus-UEL-Hyperelasticity` (BSD-3),
`jgomezc1/ABAQUS-US` (MIT), `umat_finite_viscoelasticity` (MIT),
`UMAT4COMSOL` (MIT). Even these are **not** integrated yet — they are catalogued
for planned validation. The active framework depends on **none** of them and
runs its full offline suite without any external source present.
