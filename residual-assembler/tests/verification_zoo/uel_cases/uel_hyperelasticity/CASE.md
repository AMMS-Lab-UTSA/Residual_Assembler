# CASE — Hyperelastic UEL (3D finite strain)

- **Source**: Abaqus-UEL-Hyperelasticity — https://github.com/bibekananda-datta/Abaqus-UEL-Hyperelasticity
  (mirrored: `sources/permissive/bibekanandadatta_Abaqus-UEL-Hyperelasticity`)
- **License**: BSD-3-Clause (permissive)
- **Physics**: finite-strain hyperelasticity (Neo-Hookean, Arruda-Boyce), total Lagrangian PK2/PK1
- **Element type**: 2D/3D continuum incl. **Hex8/Hex20** (full & reduced integration, F-bar)
- **User subroutine type**: UEL (standard displacement; some formulations unsymmetric)
- **Required residual backend**: `uel_direct` (Mode 3)
- **Residual class**: D — UEL / direct residual
- **Required inputs**: mesh + DOF layout + callable/exported UEL (`RHS`, `AMATRX`, `SVARS`); PROPS (mu, K, [lambda_L], nInt, matID)
- **Expected outputs**: `element_residual = -RHS`, `element_tangent = AMATRX` (may be **unsymmetric** for F-bar); global `R ~ 0`; reactions vs `RF`
- **Offline tests possible**: `uel_direct` sign/FD self-test; unsymmetric AMATRX handling check
- **Abaqus tests required**: yes — instrument/export the UEL at finite strain
- **Current status**: adapter-skeleton
- **Next missing item**: unsymmetric-tangent path + F-bar handling in `uel_direct`; confirm sign, DOF order, and state layout at large deformation
