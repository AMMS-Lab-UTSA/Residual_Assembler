# PROVENANCE — UMAT_Hyperelastic (Neo-Hookean) — REFERENCE ONLY

- **Source name**: UMAT_Hyperelastic
- **URL**: https://github.com/Sina-Taghizadeh/UMAT_Hyperelastic
- **Date reviewed**: 2026-07-13
- **License**: **GPL-3.0 (copyleft)**
- **Source category**: copyleft → **reference only**
- **Subroutine type**: UMAT
- **Physics**: Neo-Hookean finite-strain hyperelasticity (Cauchy stress + DDSDDE)
- **Element types**: continuum solids (beam-shaped test mesh)
- **Input files available**: yes (simple extension / simple shear `.inp`)
- **Abaqus validation runnable**: yes
- **Framework backend required**: `solid_c3d8_finite_strain` (Mode 2) — *if* an
  independent (non-GPL) hyperelastic material were implemented
- **Residual class**: B
- **Support status**: **reference-only** (license-limited)
- **Next missing item**: GPL-3.0 forbids copying/linking into the framework;
  re-derive Neo-Hookean stress/tangent independently from published equations if
  an in-tree hyperelastic material is ever wanted

## Isolation note
**Do NOT copy, adapt, link, or import** this code. Copyleft would relicense the
combined work. Study/benchmark reference only. Analytically validated (author
compares UMAT vs Abaqus built-in vs closed-form), so it is a good *conceptual*
check for a future, independently-written hyperelastic backend.
