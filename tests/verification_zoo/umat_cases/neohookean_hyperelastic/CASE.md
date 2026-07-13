# CASE — Neo-Hookean hyperelastic UMAT (reference only)

- **Source**: UMAT_Hyperelastic — https://github.com/Sina-Taghizadeh/UMAT_Hyperelastic
- **License**: **GPL-3.0 (copyleft)** → reference only
- **Physics**: Neo-Hookean finite-strain hyperelasticity; path-independent (elastic)
- **Element type**: continuum solid
- **User subroutine type**: UMAT
- **Required residual backend**: `solid_c3d8_finite_strain` (Mode 2) — only if an **independent** hyperelastic material is written
- **Residual class**: B — standard element + UMAT
- **Required inputs**:
  - Mode 1: mesh + exported stress field + DOF field
  - Mode 2: mesh + an independently-implemented Neo-Hookean material + PROPS (mu, K) + deformation gradient
- **Expected outputs**: residual `R`, reactions vs `RF`; the source itself is analytically validated (UMAT vs Abaqus built-in vs closed-form)
- **Offline tests possible**: Mode 1 assembly with an exported field; a **clean-room** Neo-Hookean material (re-derived from equations) could enable Mode 2 — being path-independent it does **not** need history
- **Abaqus tests required**: yes (to export fields) — but the code itself must not be used
- **Current status**: **reference-only / license-limited**
- **Next missing item**: GPL-3.0 forbids code reuse — re-implement Neo-Hookean independently if wanted. Good closed-form cross-check for a future in-tree hyperelastic backend.
