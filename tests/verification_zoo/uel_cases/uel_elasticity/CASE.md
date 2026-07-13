# CASE — Linear elastic UEL (2D/3D continuum)

- **Source**: Abaqus-UEL-Elasticity — https://github.com/bibekananda-datta/Abaqus-UEL-Elasticity
  (mirrored: `sources/permissive/bibekanandadatta_Abaqus-UEL-Elasticity`)
- **License**: BSD-3-Clause (permissive)
- **Physics**: small-strain isotropic linear elasticity
- **Element type**: Tri3/Tri6/Quad4/Quad8 (2D), Tet4/Tet10/**Hex8/Hex20** (3D)
- **User subroutine type**: UEL (standard displacement formulation)
- **Required residual backend**: `uel_direct` (Mode 3)
- **Residual class**: D — UEL / direct residual
- **Required inputs**: mesh + element DOF layout + a callable/exported UEL returning `RHS`, `AMATRX`, `SVARS`; PROPS (E, nu, nInt)
- **Expected outputs**: `element_residual = -RHS` (Abaqus sign), `element_tangent = AMATRX`; global `R ~ 0`, reactions vs `RF`; FD tangent matches AMATRX
- **Offline tests possible**: yes for the framework's `uel_direct` sign/FD self-test (already in `formulations/uel_adapter.py`); comparing to the *real* exported UEL needs instrumentation
- **Abaqus tests required**: yes — to export/instrument the UEL `RHS`/`AMATRX` (see `docs/uel_validation_plan.md`)
- **Current status**: adapter-skeleton (`uel_direct` exists; not yet wired to this specific routine)
- **Next missing item**: instrument the UEL to dump COORDS/U/PROPS/SVARS/RHS/AMATRX, then compare against `uel_direct`; confirm sign convention `RHS = -R` and DOF ordering
