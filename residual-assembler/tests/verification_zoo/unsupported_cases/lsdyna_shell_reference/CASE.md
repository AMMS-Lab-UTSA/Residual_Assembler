# CASE — LS-Dyna user shell/element (reference only, wrong solver)

- **Source**: usld/ushl LS-Dyna — https://github.com/jfriedlein/usld_LS-Dyna_Fortran , https://github.com/jfriedlein/ushl_LS-Dyna_Fortran
- **License**: unverified → reference only
- **Physics**: user element / **resultant shell** formulation basics
- **Element type**: plane-strain, axisymmetric, shell
- **User subroutine type**: LS-Dyna USLD / USHL (**not** Abaqus UEL)
- **Required residual backend**: `shell_placeholder` (contract only — not runnable)
- **Residual class**: H — reference-only due to solver/licensing
- **Required inputs**: n/a (different solver)
- **Expected outputs**: n/a
- **Offline tests possible**: no
- **Abaqus tests required**: n/a (wrong solver)
- **Current status**: **reference-only**
- **Next missing item**: an actual Abaqus shell UEL (none permissive found) or a
  native shell backend implementing the `shell_base.py` contract. Kept to document
  the shell landscape and the resultant/DOF layout of a user shell element.
