# PROVENANCE — umat_finite_viscoelasticity

- **Source name**: umat_finite_viscoelasticity
- **URL**: https://github.com/thealanjason/umat_finite_viscoelasticity
- **Date reviewed**: 2026-07-13
- **License**: MIT (permissive)
- **Source category**: permissive
- **Subroutine type**: UMAT
- **Physics**: finite viscoelasticity (1-parameter Ogden hyperelastic base + multiple relaxation mechanisms); history-dependent
- **Element types**: 3D continuum solids (e.g. C3D8)
- **Input files available**: yes (`simulation_input_files/`)
- **Abaqus validation runnable**: yes (needs Intel ifort + Abaqus to compile the UMAT)
- **Framework backend required**: `solid_c3d8_finite_strain` (Mode 2 material-replay) or `stress_driven_c3d8` (Mode 1)
- **Residual class**: B (standard element + UMAT)
- **Support status**: planned
- **Next missing item**: exported ODB stress field (Mode 1), or compiled UMAT + full increment history (Mode 2 — history-dependent, cannot final-step replay)

## Isolation note
Code is **not vendored** here. MIT permits later vendoring with attribution +
license retention; not required for any offline test. Cite Correa (2025),
Zenodo DOI 10.5281/zenodo.14801398.
