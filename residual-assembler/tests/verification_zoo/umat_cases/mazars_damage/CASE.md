# CASE — Concrete damage UMAT (Mazars)

- **Source**: Mazars_UMAT — https://github.com/marioruiarruda/Mazars_UMAT
- **License**: BSD-3-Clause **with caveat** (extra author permissions file — review before use)
- **Physics**: isotropic concrete damage (modified Mazars, fracture-energy regularized); softening / history-dependent
- **Element type**: continuum solid
- **User subroutine type**: UMAT
- **Required residual backend**: `solid_c3d8_small_strain` (Mode 2)
- **Residual class**: B — standard element + UMAT
- **Required inputs**:
  - Mode 1: mesh + exported stress field + DOF field
  - Mode 2: mesh + compiled UMAT + PROPS + increment history + damage-state STATEV
- **Expected outputs**: residual `R`; softening branch requires careful increment control; reactions vs `RF`
- **Offline tests possible**: Mode 1 assembly with an exported field (no external code run)
- **Abaqus tests required**: yes
- **Current status**: planned (**hold** — resolve the BSD-3 vs `1-Permissions_&_Instructions.txt` double-terms first)
- **Next missing item**: confirm license terms; then a single-element softening `.inp` and a Mode-1 export. Damage softening stresses the *history* requirement (Mode 2 cannot final-step replay).
