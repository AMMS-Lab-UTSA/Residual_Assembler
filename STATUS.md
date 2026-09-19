# Status

This page records what Residual_Assembler supports today and how each
capability was verified. It is for users deciding whether the program fits
their model and for reviewers checking a claim. For each workflow run end to
end with its real output, read [docs/USAGE_REPORT.md](docs/USAGE_REPORT.md);
for every quantitative claim with the command that reproduces it, read
[docs/VERIFICATION_RECORD.md](docs/VERIFICATION_RECORD.md).

One rule applies throughout: **a capability is claimed only for what was
verified against an independent reference, and the page says which
reference.** Code that runs is not, by itself, a verified result.

## Current status (2026-09-18)

### Supported workflows

| Workflow | Command | Scope | How it is verified |
|---|---|---|---|
| Sensitivities of a finished Abaqus analysis, bounded engine | `resasm request` | one homogeneous J2 C3D8/B-bar static step (NLGEOM=NO), zero fixed boundaries, ramped nodal loads, the pinned `m3_j2` provider | genuine Abaqus ODB against the analytic uniaxial J2 derivatives (relative error: E 1.7e-7, initial yield stress 8.7e-16, H 8.8e-9); the same run with every read of a Fortran source denied gives byte-identical outputs |
| History replay for any provider | `resasm history`; `resasm request` routes here automatically | small-strain C3D8 (B-bar), one static step, one user material, ramped boundaries including nonzero prescribed displacements, concentrated loads, many increments, any provider that exports `UMAT_OTI_EVAL_TOTAL` | replayed stress, state and reactions checked against the ODB at every integration point and increment; sensitivities against whole-model central differences of the ORIGINAL UMAT and against Abaqus reruns; Euler's homogeneity identity at every increment |
| Connected replay of a bounded J2 model | `resasm replay` | the pinned `m3_j2` provider on a synthetic converged FE record | whole-history derivatives against independently compiled ORIGINAL material and re-equilibrated FE paths, scaled error below `2e-6` |
| Stress-driven C3D8 assembly | `resasm assemble --mode stress-driven` | integration-point stress supplied from an export | patch tests to machine precision; unit-cube example against analytic face tractions to `1e-10` |
| Finite-strain neo-Hookean C3D8 | `resasm assemble` / `resasm sensitivity` with the bounded example configuration | stateless isotropic total neo-Hookean, first-order material parameters | independent first-Piola force (5.70e-16), complete nodal tangent FD, OTILib `du/dp` against nonlinear re-solves; reactions and reaction derivatives against one Abaqus UHYPER job |
| Direct Python residual | `resasm run` with `residual.type: python` | any order through OTILib | cubic spring against the closed form to `1e-8` at orders 1 and 2 |
| Black-box residual | `resasm run` with `residual.type: executable` | any order; the executable returns Taylor coefficients | closed form to `1e-8` at orders 1 and 2 |

Measured on the two full-size cantilevers of
[examples/cantilevers](examples/cantilevers/README.md) (evidence:
[docs/evidence/history_replay_cantilevers.md](docs/evidence/history_replay_cantilevers.md)):

| | J2 cantilever | FCC cantilever |
|---|---|---|
| mesh, DOF, integration points, increments, parameters | 1,536 C3D8, 7,497, 12,288, 40, 4 | 384 C3D8, 2,025, 3,072, 25, 10 |
| engine time, recorded state / re-equilibrated | 9.9 s / 25.2 s | 18.8 s / 71.5 s |
| worst agreement with whole-model FD of the ORIGINAL UMAT | E 6.6e-8, nu 5.8e-8, initial yield stress 7.9e-9, H 5.4e-6 | all ten parameters at most 6.7e-7 |
| homogeneity identity after re-equilibration, every increment | 1.0e-12 | 1.2e-13 |

The homogeneity identity needs no finite differences: both models are
homogeneous of degree one in their stress-dimensioned parameters, so at every
increment the parameter-weighted sensitivities of each stress and reaction sum
to the value itself, and those of each displacement sum to zero. The method
and tolerances are in [docs/REPLAY_HISTORY.md](docs/REPLAY_HISTORY.md).

Ordinary `resasm request` runs report `verified=false`: they check the replay
against the ODB but do not run the independent finite-difference reference
unless asked (`--validate`, or `resasm history --verify tangent|fd`).

### Clean installation and test suites

Both repositories were cloned from their published `main` branches into a new
directory and everything ran from those clones
([docs/evidence/final_clean_clone.md](docs/evidence/final_clean_clone.md)):

- The clean-install gate (`scripts/clean_install_gate.py`) passed all 22
  commands: provider build, `resasm request` on a genuine ODB, the
  source-denied run, both GUIs up to HTTP readiness, and the full-size J2
  cantilever (`resasm request` in 22.7 s; `resasm history --reequilibrate` in
  33.0 s, homogeneity identity 1.4e-12 of the largest term at every increment).
- The offline suite (`-m "not abaqus and not arc and not network"`) gave 494
  passed, 20 skipped and 6 failed. The six failures had one cause, an import
  check in the verification-record scripts that refused an installed
  `umat_oti` identical to the checkout's; `verification/common.py` now accepts
  a byte-identical installed copy.

[docs/COMPLETION_LEDGER.md](docs/COMPLETION_LEDGER.md) tracks every
requirement of the program and the evidence for each.

### Not supported

Refused with a named reason, not approximated:

- **Kinematics and elements:** finite-strain plasticity; element types other
  than C3D8 (C3D8R, C3D20R and C3D4 are planned). `shell_placeholder` is a
  contract only and `uel_direct` a skeleton.
- **Loads and constraints:** distributed loads (`*Dsload` is parsed but not
  applied), body forces, amplitudes, contact, and `*Equation`/MPC constraints
  (parsed but not applied).
- **Analysis structure:** several steps, several materials or instances, and
  nonzero initial state.
- **Sensitivity kinds:** load, boundary and shape sensitivities; derivatives
  above first order for finite-element models.

Other known limits:

- The Oxford crystal-plasticity UMAT in `sources/permissive/` does not compile
  under gfortran (Cray pointers and an ifort `trace()` kind mismatch). It needs
  Intel ifort and Abaqus.
- The finite-strain consistent tangent of `solid_c3d8_finite_strain` is
  analytic but approximate, pending the comparison with Abaqus `AMATRX`.
- The `resasm.yml` assembly recipe still reports `OTI-differentiate R: NO` for
  the C3D8 backends. That capability gate
  (`resasm_user/recipe.py::sensitivity_capability`) is deliberate and stays
  until that path is verified; C3D8 sensitivities come from `resasm request`,
  `resasm history` and the bounded finite-strain example instead.
- OTILib does not build natively on Windows, so the Python residual path needs
  WSL there.

## Baseline framework (July 2026)

The formulation-agnostic assembly framework and its first backend, C3D8
crystal plasticity, were established and verified offline in July 2026. That
record is kept below because its criteria still describe how the framework is
built. Where later work changed a status, the table says so.

### Framework criteria

| # | Criterion | State |
|---|---|---|
| 1 | Core makes no C3D8/CP/stress-strain/UMAT/displacement-only assumptions | Met: `core/assembler.py` only calls `Formulation.eval_element` and scatters; `dof_manager.py` builds heterogeneous per-node DOF sets |
| 2 | Crystal plasticity is only one backend | Met: CP lives in `formulations/solid_c3d8_finite_strain.py` and `materials/crystal_plasticity_adapter.py`, registered alongside the others |
| 3 | Agnosticism shown with a truss and a beam backend | Met: `truss2` (`EA/L`) and `beam2` (`PL³/3EI`); `tests/framework/test_truss2_backend.py`, `test_beam2_backend.py` |
| 4 | Shell backend contract exists | Met: `formulations/shell_base.py` and `shell_placeholder.py` (declared and documented, `supported_modes=()`, not runnable) |
| 5 | Registry inspects models and auto-selects backends | Met: `core/registry.py` and `core/diagnostics.py::inspect_model` |
| 6 | CLI and API run with minimal input | Met: `residual_core.ResidualProblem` and the `resasm` CLI assemble a real model's stress-driven residual from a field export alone |
| 7 | Requirements engine reports the minimum missing data | Met: `core/requirements.py` reports one "minimum next item" per mode |
| 8 | Docs explain adding backends | Met: `residual_core/docs/adding_a_formulation.md`, `adding_a_material.md`, `minimal_input_contract.md`, `user_interface.md` |
| 9 | Tests prove mixed-model dispatch | Met: `tests/framework/test_mixed_model_dispatch.py` (truss, beam and solid; 3/6/3 DOFs per node) |
| 10 | Existing CP verification remains intact | Met: `tests/cp_c3d8_umat/` unchanged; `test_assembler.py` reproduces the kernel bitwise in both modes, FD tangent 2.2e-13 |

Registered backends: formulations `solid_c3d8_finite_strain`,
`solid_c3d8_small_strain`, `stress_driven_c3d8`, `truss2`, `beam2`,
`nonlinear_spring1`, `nonlinear_bar1`, `shell_placeholder` (contract only) and
`uel_direct`; materials `isotropic_elastic`, `compressible_neo_hookean`, `umat`
and `crystal_plasticity`.
The generic verification ladder has Levels 0 to 7 (`core/verification.py`:
zero-field L1, rigid-body L2, FD tangent L4, reactions L5; patch tests L3 in
the C3D8 kernel).

The framework is formulation-agnostic by architecture. A formulation becomes
supported when a backend satisfying the contract is registered and verified,
not before.

### Usability hardening

| Area | Result |
|---|---|
| Minimal CLI workflow | `resasm inspect / doctor / requirements / assemble --fields / template / backends`; runnable examples under `residual_core/examples/minimal_*` |
| Requirements engine | reports only the minimum missing item and never a generic checklist; `tests/framework/test_requirements_negative.py` |
| Backend registry audit | every backend declares its required and optional inputs per mode, material interface, state requirements, tangent support and verification status; `resasm backends` and `resasm inspect --detail` render them |
| Mixed-DOF edge cases | `tests/framework/test_dof_manager_mixed.py`; the union rule is in `residual_core/docs/architecture.md` §3b |
| Public API | `tests/framework/test_public_api.py` drives everything through `residual_core.ResidualProblem`, with a static guard that no backend class is imported |
| Solver-neutral JSON | `from_neutral` and `resasm inspect model.json` round-trip nodes, elements, material and section tags, the user-material flag, boundary and load metadata, and field references (`tests/framework/test_neutral_io.py`) |

### OTILib sensitivity engine

OTILib is the hypercomplex engine behind the sensitivity paths: an external
GPLv3 dependency, built separately and never vendored. Its tests were first
run with a genuine build in WSL (Python 3.9, conda environment `pyoti`):
`test_otilib_adapter.py`, `test_otilib_spring_sensitivity.py`,
`test_otilib_fe_sensitivity.py` and `test_oti_recovery_factor.py` gave
**11 passed, 0 skipped** with `RUN_OTILIB_TESTS=1`, so a skip would have
counted as a failure. The current environment uses a Python 3.11 build
([docs/OTILIB_VENV.md](docs/OTILIB_VENV.md)).

That run includes the order-2 recovery-factor proof: for the spring
`R = k u³ − f` at `k=2, f=16, u=2`, the raw OTI coefficient for the `k²`
direction is `1/9` and the recovered derivative is `2! · 1/9 = 2/9 = d²u/dk²`,
cross-checked against a finite difference of the closed-form `u(k,f)`. Public
reports quote recovered derivatives; both conventions and
`direction_map_order<p>.json` are written to `private/`.

Reproduce with one command:

```bash
bash scripts/run_otilib_tests_wsl.sh                                   # inside WSL or Linux
powershell -ExecutionPolicy Bypass -File scripts\run_otilib_tests_wsl.ps1   # from Windows
```

The script locates conda, activates the environment holding the built
`pyoti`, exports `OTILIB_ROOT` and `PYTHONPATH`, checks that the genuine OTI
API imports, and sets `RUN_OTILIB_TESTS=1` so the tests cannot skip silently.
If anything is missing it exits non-zero and names what is missing and how to
fix it. `OTILIB_ROOT`, `OTILIB_CONDA_ENV` and `CONDA_SH` can be overridden.
`pytest` must be installed in that environment (`python -m pip install
pytest`); the script prints that command when it is absent.

On Windows without WSL the OTILib tests skip with a reason, never a silent
pass, and the Python residual path cannot run, because
`oti_global.solve_python` requires OTILib.

### Crystal-plasticity C3D8+UMAT backend

The first verified backend (`formulations/solid_c3d8_finite_strain.py` and
`materials/crystal_plasticity_adapter.py`, tests under `tests/cp_c3d8_umat/`)
was built against the criteria below. In July 2026 Abaqus was not available,
so every criterion that needs a real ODB or a compiled UMAT was built
ready-to-run and marked pending. Abaqus 2021.HF5 has since been used for the
replay engines, which close several of those gaps for UMAT-OTI providers; the
last column says which.

| # | Criterion | July 2026 | Current |
|---|---|---|---|
| 1 | The manifest classifies the downloaded examples | Verified | Unchanged |
| 2 | A permissive C3D8 CP UMAT example is parsed | Verified | Unchanged |
| 3 | Abaqus ODB fields are exported | Built, needed Abaqus | Exporters run on real ODBs: `residual_core/io/abaqus_odb_export.py` (used by `resasm request`) and `residual_core/replay/odb_export_npz.py` (used by `resasm history`) |
| 4 | Stress-driven residual reproduces Abaqus RF and free residuals | Assembly mathematics verified offline; RF comparison needed an ODB | No run of the stress-driven scripts against a real ODB is recorded; the replay engines make the equivalent check from replayed stress (row 6) |
| 5 | UMAT replay reproduces Abaqus stress and state history | Plumbing proven with a mock UMAT; real UMAT needed ifort | Replayed stress and state match the ODB at every integration point and increment for UMAT-OTI providers (J2 and FCC cantilevers); the Oxford CP UMAT itself still needs ifort |
| 6 | Residual from replayed UMAT stress matches Abaqus RF | Both halves built | Done by the replay engines for UMAT-OTI providers: reactions and the free residual are checked against the ODB at every increment |
| 7 | Tangent check passes or reports the DDSDDE mapping error | Verified, open item stated | Unchanged: the finite-strain `DDSDDE` to `AMATRX` mapping is still open |

What the July criteria established:

1. **Manifest.** All 105 `.inp` files under
   `sources/{permissive,copyleft,license-unknown}` are inventoried in
   `residual_core/example_manifest.{md,json}`. The licence invariant holds
   exactly: every copyleft or licence-unknown file is reference-only.
2. **Parsing.** `abaqus_inp_parser.py` parses `HCPnoTwin/Compression111.inp`:
   216 nodes, 125 C3D8 elements, material CPuranium (user material, 11
   constants, 125 SDV) and all boundaries, resolving the part/assembly `Set-1`
   name collision. The larger group-A jobs (DiscreteTwin, 8,035 elements;
   PyCiGen, 2,940 elements and 95 grains) also parse.
3. **Assembly mathematics, offline.** Divergence-theorem patch test on the real
   Compression111 element and a distorted hexahedron: `∫BᵀσdV` equals the
   surface-traction consistent forces to 5e-16, self-equilibrium to 1e-13. A
   linear-stress patch and global body-load test (which catches per-IP weights
   and intra-element redistribution that a uniform field cannot): 1e-13. The
   finite-strain force is frame-objective to 5e-16 and reduces exactly to
   small strain at `u = 0`.
4. **UMAT replay plumbing.** `umat_adapter_fortran/` marches a full DFGRD
   history in one process, persisting `STATEV` and the `/UMPS/` common block
   across increments. The 37-argument UMAT signature matches `umat.for`, and
   the STATEV layout matches `kmat.f`. A mock elastic UMAT shows the history
   dependence (`STATEV(35)` accumulates `0 → 1e-5 → … → 1e-4` over five
   increments; an identity increment gives zero stress).
5. **Tangent.** The small-strain element tangent `K = ∫BᵀDB` agrees with finite
   differences to 1.6e-16 (relative Frobenius). The finite-strain force
   Jacobian at fixed σ agrees with an exact analytic Jacobian to 5e-11. The
   mapping of a finite-strain UMAT `DDSDDE` to the element material tangent
   (objective rate and geometric split) that would match Abaqus `AMATRX` is
   still open; see `tests/cp_c3d8_umat/tangent_fd_check/README.md`.

Each module was cross-audited when it was built: the mathematics was
re-derived independently of the implementation and the builds re-run, and
confirmed findings were fixed and re-verified. The shared conventions are in
[residual_core/CONTRACT.md](residual_core/CONTRACT.md).

The baseline deliberately made no claim about generic elastic models, UEL or
cohesive elements (inventoried separately as group C), new crystal-plasticity
physics, or final-step parameter overloading: every replay marches the full
increment history, because `STATEV` is history-dependent.
