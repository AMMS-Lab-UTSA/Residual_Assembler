# tests/cp_c3d8_umat — verification of ONE backend

These tests verify the **C3D8 continuum + crystal-plasticity UMAT** backend
(`formulations/solid_c3d8_finite_strain.py` + `materials/crystal_plasticity_adapter.py`),
which is the *first verified formulation/material backend* of the
formulation-agnostic framework — **not** the whole framework. Framework-level
tests (that the generic assembler is formulation-agnostic) live in
[`tests/framework/`](../framework/).

Contents:

| Dir | Level(s) | Needs Abaqus? | What it checks |
|-----|----------|---------------|----------------|
| `tangent_fd_check/` | 3, 4 | no | divergence-theorem + linear-stress patch tests and the finite-difference tangent check of the C3D8 kernel (machine precision). `python tangent_fd_check/run_checks.py` |
| `stress_driven_residual/` | 5 | **yes** | procedure: export ODB fields → assemble residual externally → compare to Abaqus reactions / free residual |
| `umat_replay_vs_abaqus/` | 6, 7 | **yes** (+ifort) | procedure: replay the UMAT increment-by-increment → compare STRESS/STATEV to the ODB, then assemble a residual from replayed stress and compare to reactions |

The verification-level ladder (0–7) is defined in
[`residual_core/docs/verification_strategy.md`](../../residual_core/docs/verification_strategy.md).
This backend is the one that exercises all seven levels; other backends (built-in
material, stress-driven, UEL) use the subset that applies to them.
