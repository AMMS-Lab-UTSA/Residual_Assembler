# CASE — Cosserat / classical 2D UEL (UEL8/UEL9)

- **Source**: ABAQUS-US — https://github.com/jgomezc1/ABAQUS-US
  (mirrored: `sources/permissive/jgomezc1_ABAQUS-US`)
- **License**: MIT (permissive)
- **Physics**: 2D elasticity (classical) and **Cosserat / micropolar** (adds a micro-rotation DOF); ships classical-plasticity UMATs too
- **Element type**: 2D Quad8 (UEL8), Quad9 (UEL9)
- **User subroutine type**: UEL (+ UMATs)
- **Required residual backend**: `uel_direct` (Mode 3, 2D)
- **Residual class**: D — UEL / direct residual
- **Required inputs**: mesh + DOF layout (incl. micro-rotation for Cosserat) + callable/exported UEL (`RHS`, `AMATRX`); PROPS
- **Expected outputs**: `element_residual = -RHS`, `element_tangent = AMATRX`; global `R ~ 0`; reactions vs `RF`
- **Offline tests possible**: `uel_direct` sign/FD self-test; the extra micro-rotation DOF exercises heterogeneous DOF handling
- **Abaqus tests required**: yes — 2D UEL run + export
- **Current status**: planned (2D support + Cosserat micro-rotation DOF not yet wired)
- **Next missing item**: 2D element support in the direct adapter and a micro-rotation DOF type; then compare `RHS`/`AMATRX`
