# CASE — C3D8 crystal-plasticity, stress-driven (internal, verified)

- **Source**: internal — `sources/permissive/ngrilli_Oxford_Crystal_Plasticity` (MIT) via the framework's own `stress_driven_c3d8` backend
- **License**: MIT (permissive)
- **Physics**: crystal plasticity (one example backend), assembled **stress-driven** (no material update)
- **Element type**: C3D8 (8-node hex)
- **User subroutine type**: UMAT (Grilli CP) — but here the residual is built from an **exported** stress field, so the UMAT is not called
- **Required residual backend**: `stress_driven_c3d8` (Mode 1)
- **Residual class**: F — stress-driven verification candidate
- **Required inputs**: mesh + exported integration-point Cauchy stress `S` (Voigt) + DOF field
- **Expected outputs**: `f_int = sum_k B^T sigma_k dV`; free-DOF `R ~ 0`; reactions vs `RF`
- **Offline tests possible**: **yes — offline-verified today.** `tests/framework/test_assembler.py` reproduces the verified C3D8 kernel bitwise from a stress field; `residual_core/examples/minimal_c3d8_stress_driven` assembles a nonzero residual from a supplied field with the `resasm` CLI
- **Abaqus tests required**: only to compare against a *real* Abaqus ODB export on a spatially-varying state (Level 5); the assembly math itself is verified offline
- **Current status**: **implemented/offline-verified** (Abaqus ODB cross-check pending)
- **Next missing item**: a real ODB stress export on a non-uniform state to confirm Abaqus' integration-point ordering (see `docs/limitations.md`)

## Why this is the template for every external UMAT
Any external UMAT case (viscoelastic, J2, damage, hyperelastic) reduces to **this
same offline-verifiable assembly** once its stress field is exported — the
material model becomes irrelevant to the residual assembly. That is the
formulation-agnostic claim in action.
