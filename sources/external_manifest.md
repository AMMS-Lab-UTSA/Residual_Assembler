# External Verification Candidates — Manifest

Review date: **2026-07-13**. Machine-readable twin: [external_manifest.json](external_manifest.json).

This is a **map of the landscape**, not a claim of support. It catalogues public
Abaqus user-subroutine (UMAT/UEL) examples across diverse physics so the
framework's *formulation-agnostic* claim can be checked against many external
material models and user elements — not only crystal plasticity.

## Acquisition & isolation policy (read first)

- **Code is NOT vendored** into this repository. Only provenance metadata (name,
  URL, license, physics, element types, backend) is recorded. This is deliberate:
  it keeps the permissive core clean and avoids copying restricted or
  unknown-license code.
- **Copyleft (GPL/AGPL)**, **custom/unknown-license**, and **different-solver**
  sources are **reference-only**: never copied, adapted, linked, or imported.
  Recorded under [reference_only/](reference_only/).
- **Permissive (MIT / BSD-3 / Apache)** sources *may* be vendored later **with
  attribution and license retention**. Candidates not yet adopted are recorded
  under [external_candidates/](external_candidates/) and are not required for any
  offline test. One already-adopted permissive submodule *is* required:
  `permissive/ngrilli_Oxford_Crystal_Plasticity` supplies the C3D8 mesh used by
  the assembler/recipe/neutral-IO tests. Fetch it with
  `./scripts/init_permissive_sources.sh`; see [SUBMODULES.md](SUBMODULES.md).
- No license = all rights reserved = **reference-only** by default.

> "On GitHub" does not mean "usable". License first, always.

## Diversity (acceptance criterion 4)

Physics covered: linear elasticity, finite-strain hyperelasticity, **J2
plasticity**, viscoelasticity, damage, cohesive/traction-separation, and
**coupled chemo-mechanics** — plus crystal plasticity, which remains **one
backend among many**, not the center.

## UMAT candidates (≥5)

| # | Name | License | Tier | Physics | Element | .inp | Backend / Mode | Status |
|---|------|---------|------|---------|---------|------|----------------|--------|
| U1 | [umat_finite_viscoelasticity](https://github.com/thealanjason/umat_finite_viscoelasticity) | **MIT** | permissive | viscoelasticity (Ogden+relaxation) | C3D8 solids | yes | `solid_c3d8_finite_strain` · Mode 2 / Mode 1 | planned |
| U2 | [UMAT4COMSOL](https://github.com/sergiolucarini/UMAT4COMSOL) (elastoplastic) | **MIT** | permissive | **J2 von-Mises plasticity** (+ neo-Hookean, CP) | solids | no (COMSOL) | `solid_c3d8_*` · Mode 2 | planned |
| U3 | [Mazars_UMAT](https://github.com/marioruiarruda/Mazars_UMAT) | **BSD-3** ⚠ | permissive | concrete **damage** (Mazars) | solids | partial | `solid_c3d8_small_strain` · Mode 2 | planned (see ⚠) |
| U4 | [UMAT_Hyperelastic](https://github.com/Sina-Taghizadeh/UMAT_Hyperelastic) | **GPL-3.0** | copyleft | hyperelasticity (Neo-Hookean) | solids | yes | `solid_c3d8_finite_strain` · Mode 1/2 | reference-only |
| U5 | [Bilinear-CZM-UMAT](https://github.com/harshaa765/Bilinear-CZM-UMAT) | **GPL-3.0** | copyleft | **cohesive** (bilinear TSL) | 2D cohesive | yes | none yet · Mode 2 | unsupported / reference-only |

⚠ **U3 (Mazars)**: LICENSE is BSD-3, but the repo also ships a
`1-Permissions_&_Instructions.txt` with author-imposed conditions. Review those
terms before any reuse; treated as permissive-with-caveat.

## UEL candidates (≥5)

| # | Name | License | Tier | Physics | Element | .inp | Backend / Mode | Status |
|---|------|---------|------|---------|---------|------|----------------|--------|
| E1 | [Abaqus-UEL-Elasticity](https://github.com/bibekananda-datta/Abaqus-UEL-Elasticity) | **BSD-3** | permissive | linear elasticity | 2D/3D Tri/Quad/Tet/**Hex8/20** | yes | `uel_direct` · Mode 3 | adapter-skeleton |
| E2 | [Abaqus-UEL-Hyperelasticity](https://github.com/bibekananda-datta/Abaqus-UEL-Hyperelasticity) | **BSD-3** | permissive | finite-strain hyperelastic (PK2/PK1) | 2D/3D incl **Hex8/20** | yes | `uel_direct` · Mode 3 | adapter-skeleton |
| E3 | [Abaqus-UEL-Hydrogel](https://github.com/bibekananda-datta/Abaqus-UEL-Hydrogel) | **custom** (non-SPDX) | unknown | **coupled chemo-mechanics** | coupled (u+μ) | yes | coupled-field · Mode 5 (n/impl) | reference-only |
| E4 | [ABAQUS-US](https://github.com/jgomezc1/ABAQUS-US) (UEL8/9 + Cosserat) | **MIT** | permissive | elasticity, **Cosserat**, plasticity UMAT | 2D Quad8/9 | yes | `uel_direct` · Mode 3 (2D) | planned |
| E5 | [Abaqus-UEL-Subroutine](https://github.com/CAEAssistant-Group/Abaqus-UEL-Subroutine) | **MIT** (demo) | reference-only | generic UEL skeleton | unspecified | no | `uel_direct` (concept) | reference-only (incomplete) |
| E6 | [usld/ushl LS-Dyna](https://github.com/jfriedlein/usld_LS-Dyna_Fortran) | unverified | reference-only | user element / **shell** basics | shell, plane, axi | no | n/a (LS-Dyna) | reference-only (wrong solver) |

## Documented gaps (honest)

- **Truss/bar UEL (Abaqus)** — none permissive found; covered by native `truss2`
  (Mode 4).
- **Beam/frame UEL (Abaqus)** — none permissive found; covered by native `beam2`
  (Mode 4).
- **Standalone thermal/diffusion UEL (Abaqus)** — closest is the coupled hydrogel
  UEL (E3); Mode-5 coupled DOFs are designed-for but not implemented.

## Already-mirrored permissive sources (for reference)

E1, E2, E4 are already present under [permissive/](permissive/) from the earlier
inventory ([inventory.md](inventory.md)). They are listed here again in the
external candidate framing so the verification-zoo cards can point at one place.

See [LICENSE_AUDIT.md](LICENSE_AUDIT.md) for the per-source license determination
and the integration decision for each.
