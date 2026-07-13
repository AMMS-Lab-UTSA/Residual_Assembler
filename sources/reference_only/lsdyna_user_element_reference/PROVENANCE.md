# PROVENANCE — usld/ushl LS-Dyna user element — REFERENCE ONLY

- **Source name**: usld_LS-Dyna_Fortran / ushl_LS-Dyna_Fortran
- **URL**: https://github.com/jfriedlein/usld_LS-Dyna_Fortran ,
  https://github.com/jfriedlein/ushl_LS-Dyna_Fortran
- **Date reviewed**: 2026-07-13
- **License**: unverified
- **Source category**: **reference-only (different solver)**
- **Subroutine type**: LS-Dyna user element (USLD) / user shell (USHL) — **not Abaqus**
- **Physics**: user-defined element / resultant shell formulation basics
- **Element types**: plane-strain, axisymmetric, **shell**
- **Input files available**: no Abaqus `.inp` (LS-Dyna keyword decks)
- **Abaqus validation runnable**: no (wrong solver)
- **Framework backend required**: n/a
- **Residual class**: H (reference-only due to solver/licensing)
- **Support status**: **reference-only**
- **Next missing item**: LS-Dyna ≠ Abaqus; conceptual reference for the
  structure of a user shell/element residual and its resultant/DOF layout

## Isolation note
Catalogued to document the **shell** and user-element landscape (the framework's
`shell_placeholder` contract has no runnable external Abaqus example yet). Not an
Abaqus source; do not integrate.
