# CASE — J2 (von-Mises) plasticity UMAT

- **Source**: UMAT4COMSOL (`elastoplastic/`) — https://github.com/sergiolucarini/UMAT4COMSOL
- **License**: MIT (permissive)
- **Physics**: rate-independent J2 von-Mises elastoplasticity; history-dependent (plastic strain / hardening state)
- **Element type**: continuum solid (map to C3D8)
- **User subroutine type**: UMAT (standard Abaqus interface)
- **Required residual backend**: `solid_c3d8_small_strain` or `solid_c3d8_finite_strain` (Mode 2)
- **Residual class**: B — standard element + UMAT
- **Required inputs**:
  - Mode 1: mesh + exported stress field + DOF field
  - Mode 2: mesh + compiled UMAT + PROPS (E, nu, sigma_y, hardening) + increment history + STATEV
- **Expected outputs**: residual `R`, free-DOF `R ~ 0`, reactions vs `RF`; in Mode 2 also STRESS/STATEV vs ODB
- **Offline tests possible**: Mode 1 assembly if stress export provided; a Python J2 return-map could be **independently** written for Mode 2 (do not copy the source)
- **Abaqus tests required**: yes (build the `.inp` — repo ships no Abaqus deck — and run/export)
- **Current status**: planned
- **Next missing item**: author a single-element Abaqus `.inp` for this UMAT; then compiled UMAT (Mode 2) or exported stress (Mode 1). This is the key **non-crystal-plasticity** plasticity example.
