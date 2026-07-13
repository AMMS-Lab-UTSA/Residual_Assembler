# PROVENANCE — UMAT4COMSOL (J2 plasticity)

- **Source name**: UMAT4COMSOL
- **URL**: https://github.com/sergiolucarini/UMAT4COMSOL
- **Date reviewed**: 2026-07-13
- **License**: MIT (permissive)
- **Source category**: permissive
- **Subroutine type**: UMAT (standard Abaqus format) + a C/COMSOL wrapper
- **Physics**: J2 von-Mises elastoplasticity (`elastoplastic/`), neo-Hookean
  hyperelasticity (`neohookean/`), crystal plasticity (`crystalplasticity/`),
  coupled (`coupled/`). The **J2 plasticity** routine is the diverse,
  non-crystal-plasticity example of interest.
- **Element types**: continuum solids
- **Input files available**: no Abaqus `.inp` (repo is COMSOL-oriented; ships
  `result_files/`). A single-element Abaqus `.inp` must be authored to validate.
- **Abaqus validation runnable**: yes for the UMAT itself (standard interface)
- **Framework backend required**: `solid_c3d8_small_strain` / `solid_c3d8_finite_strain` (Mode 2)
- **Residual class**: B
- **Support status**: planned
- **Next missing item**: author an Abaqus single-element `.inp` for the J2 UMAT;
  then compiled UMAT (Mode 2) or exported stress (Mode 1)

## Isolation note
Code is **not vendored**. MIT; may be vendored later with attribution. Cite
Lucarini & Martinez-Paneda, *Advances in Engineering Software* (2024).
