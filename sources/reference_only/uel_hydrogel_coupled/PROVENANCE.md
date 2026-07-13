# PROVENANCE — Abaqus-UEL-Hydrogel (coupled) — REFERENCE ONLY

- **Source name**: Abaqus-UEL-Hydrogel
- **URL**: https://github.com/bibekananda-datta/Abaqus-UEL-Hydrogel
- **Date reviewed**: 2026-07-13
- **License**: **custom / non-SPDX** ("View license" — bespoke `LICENSE.md`)
- **Source category**: unknown/custom → **reference only** (until terms confirmed)
- **Subroutine type**: UEL
- **Physics**: coupled chemo-mechanics of hydrogels (species diffusion coupled to
  finite deformation) — a genuine **multi-field / coupled-DOF** example
- **Element types**: coupled continuum (displacement + chemical-potential DOFs)
- **Input files available**: yes (`test/`)
- **Abaqus validation runnable**: yes
- **Framework backend required**: **coupled-field formulation (Mode 5)** — not
  implemented (the DOF manager supports heterogeneous DOFs, but no coupled
  constitutive/DOF wiring exists yet)
- **Residual class**: E (coupled-field formulation)
- **Support status**: **reference-only**
- **Next missing item**: confirm the custom LICENSE terms permit reuse; implement
  coupled-field DOFs (u + μ) and the coupled weak form

## Isolation note
License is a **bespoke file**, not a recognized SPDX identifier → treat as
all-rights-reserved until reviewed. This is the best available coupled/multi-field
UEL example and a strong motivator for the Mode-5 roadmap item. Author's simpler
BSD-3 UELs (Elasticity, Hyperelasticity) are the permissive starting points
instead.
