# CASE — Finite viscoelasticity UMAT

- **Source**: umat_finite_viscoelasticity — https://github.com/thealanjason/umat_finite_viscoelasticity
- **License**: MIT (permissive)
- **Physics**: finite viscoelasticity (Ogden hyperelastic base + relaxation); history-dependent
- **Element type**: 3D continuum solid (C3D8)
- **User subroutine type**: UMAT
- **Required residual backend**: `solid_c3d8_finite_strain` (Mode 2) or `stress_driven_c3d8` (Mode 1)
- **Residual class**: B — standard element + UMAT
- **Required inputs**:
  - Mode 1: mesh + exported integration-point stress field (ODB) + DOF field
  - Mode 2: mesh + compiled UMAT + PROPS + **full increment history** + dtime per step + initial STATEV
- **Expected outputs**: global residual `R`; free-DOF `R ~ 0` at equilibrium; reactions matching Abaqus `RF`
- **Offline tests possible**: Mode 1 residual assembly *if* a stress export is provided (no Abaqus needed to assemble). No external code executed.
- **Abaqus tests required**: yes — to produce the ODB stress export (Mode 1) or to compile/run the UMAT (Mode 2, needs Intel ifort + Abaqus)
- **Current status**: planned
- **Next missing item**: an exported ODB stress field (Mode 1) — then this is offline-verifiable; Mode 2 additionally needs the compiled UMAT and the increment sequence (history-dependent: no final-step-only replay)
