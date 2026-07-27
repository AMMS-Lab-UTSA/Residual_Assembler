# Verification Strategy — the generalized Level 0-7 ladder

Verification is generalized into a ladder that any formulation backend inherits
through the formulation-agnostic core. The levels are defined in
`core/verification.py` (module docstring) and each backend lists which apply in
its `verification_levels` attribute. The governing rule (from `STATUS.md`): **do
not claim success because the code runs** — everything that does not need Abaqus
is verified here to machine precision; everything that needs Abaqus is built
ready-to-run and labelled pending, never passing.

The ladder:

```
0  parser / model inventory
1  zero-field / zero-load residual
2  rigid-body / null-mode
3  manufactured solution / patch test
4  finite-difference tangent
5  solver comparison vs Abaqus (exported fields -> external residual -> reactions/free residuals)
6  material replay (standalone material update -> compare stress/state history to ODB)
7  full residual replay (material replay + assembly -> compare reactions)
```

---

## Level 0 — parser / model inventory
- **Proves:** the model was ingested correctly (right node/element counts,
  element types, material user-flag / constants / SDV, BC sets) so everything
  above it operates on the true mesh.
- **Needs Abaqus?** No.
- **Realized in code:** `io/abaqus_inp_parser.parse_inp` (self-test asserts
  Compression111 = 216 nodes, 125 C3D8, material CPURANIUM user_material,
  11 constants, 125 depvar); the example manifest classification
  (`residual_core/example_manifest.md`).

## Level 1 — zero-field / zero-load residual
- **Proves:** at `U = 0` with a zero stress field (or a material returning zero at
  zero strain), the assembled internal force vanishes, so `R = -F_external`; no
  spurious forces from the assembly/scatter.
- **Needs Abaqus?** No.
- **Realized in code:** `core/verification.zero_field_residual(assembler)` — runs
  `assembler.assemble(zeros)` and returns `max|R|`. Exercised by
  `tests/framework/test_assembler.py::test_C_zero_field`.

## Level 2 — rigid-body / null-mode
- **Proves:** rigid-body motion produces no strain and hence no internal force
  (translation = a null mode of `R`); under finite strain, a superposed rigid
  rotation must leave the material response frame-objective.
- **Needs Abaqus?** No.
- **Realized in code:** `core/verification.py::rigid_body_translation` is the
  dedicated Level-2 routine — a uniform rigid-body translation must give zero
  internal force for material-driven backends (verified through the framework:
  `max|R| ~ 6e-11` on the elastic Compression111 model). What is NOT yet in an
  automated test is a rigid-ROTATION objectivity check: the finite-strain kernel's
  frame-objectivity (force rotates with R) and its exact reduction to small strain
  at `u = 0` are asserted in `residual_core/README.md` / `STATUS.md` (reported
  5e-16, and re-derived by audit) but are not yet exercised by a standing test.
  Level 2 is declared applicable to the solid backends but is the least
  independently instrumented rung.

## Level 3 — manufactured solution / patch test
- **Proves:** the element reproduces exact known nodal forces for a prescribed
  (constant or linear) stress field — validating the B-matrix, `detJ`,
  quadrature weights, node ordering, and (for linear fields) per-IP weights and
  the sigma<->point pairing.
- **Needs Abaqus?** No.
- **Realized in code:** `formulations/c3d8_kernel.py` +
  `tests/cp_c3d8_umat/tangent_fd_check/run_checks.py`:
  - divergence-theorem patch test (`surface_traction_nodal_forces` vs
    `element_internal_force_small_strain`), real Compression111 element + a
    distorted hex, random + uniaxial stress — measured `~5e-16`;
  - linear-stress patch (`linear_stress_target`, per-IP weights + pairing) —
    `~4e-14`;
  - uniform uniaxial sanity on a unit cube — `~1e-14`.
  Does NOT probe the Abaqus IP *ordering* (offline stresses are sampled at our own
  points).

## Level 4 — finite-difference tangent
- **Proves:** the assembled tangent `K` equals `dR/dU` — catches
  DDSDDE->element-tangent mapping errors for any formulation/material.
- **Needs Abaqus?** No.
- **Realized in code:** `core/verification.finite_difference_tangent(assembler, U,
  ...)` — central-differences the assembled residual w.r.t. `U` (optionally a DOF
  subset) and compares to the assembled `K`; returns `(abs_err, rel_fro,
  max_entry, K, Kfd)`. Backend-independent (runs through `Assembler`). Also, the
  kernel-level checks `check3a_linear_tangent` (small-strain `K = sum B0^T D B0`,
  relative Frobenius `~1.6e-16`) and `check3b_finite_strain` (finite-strain force
  FD vs the exact fixed-sigma Jacobian `force_tangent_fixed_sigma`, `~5e-11`) in
  `run_checks.py`. The finite-strain *material* tangent (DDSDDE mapping) is the
  deferred item, documented in
  `tests/cp_c3d8_umat/tangent_fd_check/README.md`.

## Level 5 — solver comparison vs Abaqus
- **Proves:** the residual assembly reproduces Abaqus equilibrium — feed Abaqus'
  exported integration-point stress `S` and displacement `U` through the external
  assembler; free-DOF residual must be ~0 and the prescribed-DOF internal force
  must equal the Abaqus reaction `RF` up to one global sign.
- **Needs Abaqus?** **Yes** (a real ODB).
- **Realized in code:** the Mode-1 driver
  `residual_core/stress_driven_residual.py` + `core/verification.reaction_check`
  (splits `R` into free residual and reactions, compares to a reference `RF` with
  sign). Fields come from `io/abaqus_odb_export.py` (ODB -> `fields.json`).
  Procedure and pass criteria: `tests/cp_c3d8_umat/stress_driven_residual/README.md`.

## Level 6 — material replay vs Abaqus ODB
- **Proves:** the standalone material update reproduces the Abaqus stress/state
  history — march the deformation-gradient history increment-by-increment and
  compare replayed `STRESS` to ODB `S` and `STATEV` to ODB `SDV`.
- **Needs Abaqus?** **Yes** (ODB), and for the real CP UMAT also **Intel ifort +
  Abaqus/MKL** to compile it.
- **Realized in code:** the standalone replay orchestrator
  `residual_core/umat_adapter_fortran/umat_replay.py` (drives the compiled Fortran
  `umat_driver`, persisting `STATEV` + the `/UMPS/` common block across
  increments; never a final-step jump). Procedure:
  `tests/cp_c3d8_umat/umat_replay_vs_abaqus/README.md`. Proven end-to-end today
  only with the **mock** elastic UMAT (`umat_replay.py --dry-run`); the real CP
  object is blocked on the toolchain.

## Level 7 — full residual replay
- **Proves:** the composition of Level 6 into Level 5 — assemble a residual from
  the *replayed* stress field and compare to `RF`, so the material replay and the
  assembler simultaneously agree with Abaqus.
- **Needs Abaqus?** **Yes** (ODB + ifort for the real UMAT).
- **Realized in code:** collect replayed per-IP stress into
  `sigma_all[eid] -> (8,6)` and reuse the Mode-1 assembler
  (`c3d8_kernel.assemble_global_internal_force(..., mode='finite')`) with the
  Level-5 pass criteria. Both halves are built; the end-to-end number awaits
  Abaqus + ifort.

---

## Per-backend applicability and current status

Levels marked applicable come from each backend's `verification_levels`. Status is
honest against `STATUS.md`: offline levels (0-4) are verified to machine precision;
Abaqus-dependent levels (5-7) are built and ready but not executed here (Abaqus is
not installed; `which abaqus` -> none). Level 2 is applicable-but-not-instrumented
as noted above.

| Level | solid_c3d8_finite_strain | solid_c3d8_small_strain | stress_driven_c3d8 | uel_direct |
|-------|--------------------------|-------------------------|--------------------|-----------|
| 0 parser/inventory        | applies — passing offline | applies — passing offline | applies — passing offline | n/a |
| 1 zero-field residual     | applies — passing offline | applies — passing offline | applies — passing offline | applies — passing (self-test) |
| 2 rigid-body/null-mode    | applies — not instrumented (see above) | applies — not instrumented | n/a | n/a |
| 3 patch/manufactured      | applies — passing offline (5e-16 / 4e-14) | applies — passing offline | applies — passing offline | n/a |
| 4 FD tangent              | applies — force FD passing; material DDSDDE mapping deferred | applies — passing offline (rel 1.6e-16) | n/a (no tangent in Mode 1) | applies — passing (self-test FD) |
| 5 solver vs Abaqus (RF)   | applies — built, awaits Abaqus run | applies — built, awaits Abaqus run | applies — built, awaits Abaqus run | applies — built, awaits a UEL + run |
| 6 material replay vs ODB  | applies — plumbing proven (mock); real UMAT awaits ifort+Abaqus | n/a | n/a | n/a (no material to replay) |
| 7 full residual replay    | applies — both halves built; awaits 5+6 | n/a | n/a | applies (assembly half) — awaits a run |

Notes:
- "n/a" = the level is not in that backend's declared `verification_levels`.
  `stress_driven_c3d8` declares `(0,1,3,5)` (no tangent, no material replay);
  `uel_direct` declares `(1,4,5,7)` (a direct-residual element has no IP field to
  replay, so 6 and the material half of 7 do not apply);
  `solid_c3d8_small_strain` declares `(0,1,2,3,4,5)` (small-strain elastic, no
  Abaqus CP replay path).
- The Abaqus **C3D8 IP ordering** (`ABAQUS_C3D8_GAUSS`) is high-confidence
  (CalculiX-corroborated lexicographic) but is provably not testable offline; it
  must be confirmed against the first ODB with a spatially-varying stress before
  trusting Level 5 on a non-uniform state.
