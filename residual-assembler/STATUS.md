# STATUS — Model-Agnostic Residual Assembly Framework (CP is one example backend)

Date: 2026-07-10. Honest accounting. The rule "**do not claim success because the
code runs**" is applied throughout: each claim says exactly what was verified and
how, and what is still pending a real Abaqus run. Two parts:

- **Part I — the refactor** into a formulation-agnostic framework (below).
- **Part II — the crystal-plasticity C3D8+UMAT backend** (the original Step-10
  criteria), unchanged and still verified, now living as **one example backend**
  behind the interfaces.

Offline suite green: framework (assembler + truss2 + beam2 + mixed-model
dispatch), both C3D8 solid backends, UEL adapter, UMAT/CP material adapters via a
mock UMAT, parser, ODB-export compile, neutral-model IO round-trip, runnable API
examples, manifest.

> **Note (Part 0 supersedes the "Abaqus not installed" caveat for THIS machine).**
> The dev/build environment where Parts I–II were written had no Abaqus. The user's
> workstation does — see Part 0. The ODB-comparison beats are no longer all ⏳: the
> elastic stress-driven reaction match and an Abaqus-verified sensitivity run are now
> ✅ (dated 2026-07-15).

---

## Part 0 — Live Abaqus verification on this workstation (2026-07-15)

Abaqus 2024 is installed here (`C:\SIMULIA\Commands\abaqus.bat`) with a working
FlexNet academic license. Two beats that Parts I–II left ⏳ for lack of Abaqus are now
executed and passing for the elastic case; a full-picture write-up is in
[METHODOLOGY.md](METHODOLOGY.md).

**Assembly of R vs live Abaqus reactions — ✅ (elastic).** A single elastic C3D8,
1/8-symmetry uniaxial `nlgeom=NO` job solved in Abaqus, exported with
`scripts/extract_odb_fields.py`, re-assembled here:

- `||R_free||` = 1.24e-14 (free-DOF equilibrium)
- `rel |R - RF|` = **7.82e-17** (assembled reaction vs Abaqus `RF`, sign `R = +RF`)

Driven by a **single job file** with the new `resasm verify-job` command; runnable
example: `examples/abaqus_elastic_c3d8/`.

**Sensitivity method vs Abaqus finite differences — ✅ (elastic).** The residual
method (`R^(1)` -> `K U^(1) = -R^(1)` -> chain rule) produced `du/dE` and `dσ11/dE`
matching Abaqus (solved at `E±`, central differenced) to ~1e-5 — the FD truncation,
the engine being exact. Example: `examples/residual_sensitivity_c3d8/`.

**Still ⏳ against Abaqus:** non-uniform stress (C3D8 IP ordering), finite-strain
(`nlgeom=YES`), and the crystal-plasticity UMAT. The UMAT's `d sigma / da` comes from
the **OTIS source-transformation** tool; the connection and its pending items (a
parameter-seed transform, and a local Windows build of the OTIS Fortran — Abaqus here
wants classic `ifort`, only `ifx` is present) are documented in
[docs/otis_umat_connection.md](docs/otis_umat_connection.md).

---

## Part I — Refactor success criteria

| # | Criterion | State |
|---|-----------|-------|
| 1 | Core makes no C3D8/CP/stress-strain/UMAT/displacement-only assumptions | ✅ `core/assembler.py` only calls `Formulation.eval_element` + scatters (grep-clean); `dof_manager.py` builds heterogeneous per-node DOF sets |
| 2 | Crystal plasticity is only one backend | ✅ CP lives in `formulations/solid_c3d8_finite_strain.py` + `materials/crystal_plasticity_adapter.py`, registered alongside others |
| 3 | Proven agnostic with a truss AND a beam backend | ✅ `truss2` (`EA/L`), `beam2` (`PL³/3EI`) verified; `tests/framework/test_truss2_backend.py`, `test_beam2_backend.py` |
| 4 | Shell backend contract exists | ✅ `formulations/shell_base.py` + `shell_placeholder.py` (declared, documented, `supported_modes=()` = not runnable) |
| 5 | Registry inspects models & auto-selects backends | ✅ `core/registry.py` + `core/diagnostics.py::inspect_model` (per-element auto-selection, reachable modes) |
| 6 | CLI + API run with minimal input | ✅ `residual_core.ResidualProblem` + `resasm` CLI (`ui/`); assemble a real model's stress-driven residual with only a field export |
| 7 | Requirements engine reports minimum missing data | ✅ `core/requirements.py` (per-mode minimum inputs, single "minimum next item") |
| 8 | Docs explain adding backends | ✅ `docs/adding_a_formulation.md`, `docs/adding_a_material.md`, `docs/minimal_input_contract.md`, `docs/user_interface.md` |
| 9 | Tests prove mixed-model dispatch | ✅ `tests/framework/test_mixed_model_dispatch.py` (truss + beam + solid, 3/6/3 DOFs/node, physics-blind scatter) |
| 10 | Existing CP verification remains intact | ✅ `tests/cp_c3d8_umat/` unchanged, identical numbers; `test_assembler.py` reproduces the kernel **bitwise** (Mode 1 + Mode 2), FD tangent 2.2e-13 |
| 11 | README no longer CP-centered | ✅ `residual_core/README.md` reframed as a model-agnostic framework; CP listed as one example backend |

Backends registered: formulations `solid_c3d8_finite_strain`,
`solid_c3d8_small_strain`, `stress_driven_c3d8`, `truss2`, `beam2`,
`shell_placeholder` (contract only), `uel_direct`; materials `isotropic_elastic`,
`umat`, `crystal_plasticity`. Verification ladder generalized to Levels 0–7
(`core/verification.py`: zero-field L1, rigid-body L2, FD-tangent L4, reactions L5;
patch tests L3 in the CP kernel).

**Honest claim.** The framework is formulation-agnostic *by architecture*. Each
formulation becomes supported when a backend satisfying the formulation contract
is registered and verified.

What the refactor deliberately did **not** do: no OTI/HYPAD, no 150-parameter
model, no non-mechanical coupled-field DOFs (thermal/pressure), no claim of
arbitrary-formulation support.

Continuum-solid coverage is **six formulations** — C3D8, C3D8R, C3D20, C3D20R,
C3D4, C3D10 — all validated for **assembly** (shape-function patch test and
`rel(R−RF)` reaction equilibrium ≤ 4.7e-09; see
[results/figures/element_validation.json](results/figures/element_validation.json)).
Of these, **four have full-solve verification** (the assembled `d(σ_vM)/dp`
against finite differences of a full perturbed solve): C3D8 (≤2.7e-08),
C3D20 (≤8.6e-08), C3D4 (≤1.0e-08), C3D10 (≤5.7e-08). **C3D8R and C3D20R are
assembly-only** (`full_solve_rel: null`) — their reduced integration is exercised
in assembly but not yet closed with an independent full-solve reference. (An
earlier revision of this note said "C3D8 only, C3D8R/C3D20R/C3D4 planned"; that
predates the element sweep and is superseded by the numbers above.) Truss/beam
are small-strain linear proof backends, not production elements; the shell
backend is a contract, not runnable.

### Part I-b — User-facing agnosticism audit & interface hardening

Hardening pass focused on usability by a non-developer (no low-level object
instantiation), verified offline:

| Area | Result |
|---|---|
| Minimal CLI workflow | ✅ `resasm inspect / doctor / requirements / assemble --fields / template / backends`; four runnable examples under `residual_core/examples/minimal_*` (model + command + expected output + auto-detected vs. manual notes) |
| Requirements engine | ✅ reports only the **minimum** missing item, never crashes, never a generic checklist; `tests/framework/test_requirements_negative.py` (beam/no-section, solid/no-field, UMAT/no-source, shell-placeholder, UEL/no-callable) |
| Backend registry audit | ✅ every backend declares `required_inputs_by_mode`, `optional_inputs_by_mode`, `material_interface_needed`, `state_requirements`, `tangent_support`, `verification_status`; `resasm backends` + `resasm inspect --detail` render per-element selection/status/limitations/modes/next-step |
| Mixed-DOF edge cases | ✅ `tests/framework/test_dof_manager_mixed.py` (truss-only, beam-only, truss∪beam, beam∪solid); union rule documented in `docs/architecture.md` §3b |
| Public API | ✅ `tests/framework/test_public_api.py` drives everything through `residual_core.ResidualProblem` with a static guard that no backend/internal classes are imported |
| Solver-neutral JSON | ✅ `from_neutral` + `resasm inspect model.json`; round-trips nodes, elements/types, material/section tags, `user_material` flag, boundary/load metadata, and optional `field_refs` (`tests/framework/test_neutral_io.py`) |
| Docs / README | ✅ README leads with "formulation-agnostic residual assembly framework"; CP listed under example backends only |

---

### Part I-c — OTILib (hypercomplex sensitivity engine): PASS PATH VERIFIED IN WSL

**Status: ✅ OTILib pass path verified in WSL.** This is a real numeric run, not a skip.

| item | value |
|---|---|
| OTILib build | `/root/otilib` (source build, GPLv3, external — **not** vendored, **not** a dependency) |
| conda env | `pyoti` at `/root/miniconda3/envs/pyoti` (Python 3.9.23, numpy 2.0.2) |
| built extension | `/root/otilib/build/pyoti/sparse.cpython-39-x86_64-linux-gnu.so` |
| adapter detection | `otilib_available()` → `True`; `api_module` = `pyoti.sparse` (genuine OTI API verified, not the unrelated PyPI squat) |
| tests run | `tests/framework/test_otilib_adapter.py`, `test_otilib_spring_sensitivity.py`, `test_otilib_fe_sensitivity.py`, `test_oti_recovery_factor.py` |
| result | **11 passed** (0 skipped), with `RUN_OTILIB_TESTS=1` so the skip guard is disabled and a skip would have been a failure |

Included in that run: the **order-2 recovery-factor proof** — for the spring
`R = k u³ − f` at `k=2, f=16, u=2`, the raw OTI coefficient for the `k²` direction
is `1/9` and the recovered derivative is `2! · 1/9 = 2/9 = d²u/dk²`, cross-checked
against a finite difference of the closed-form `u(k,f)`. Public reports quote the
**recovered derivatives**; both conventions plus `direction_map_order<p>.json` ship
in `private/`.

Reproduce with **one command**:

```bash
# from inside WSL / Linux
bash scripts/run_otilib_tests_wsl.sh
```
```powershell
# from Windows (launches WSL for you)
powershell -ExecutionPolicy Bypass -File scripts\run_otilib_tests_wsl.ps1
```

`scripts/run_otilib_tests_wsl.sh` makes activation **deterministic**: it locates
conda, activates the env holding the built `pyoti`, exports `OTILIB_ROOT` (+
`PYTHONPATH`), asserts the genuine OTI API actually imports, and exports
`RUN_OTILIB_TESTS=1` so the tests **cannot silently skip**. If any piece is
missing it exits non-zero and names exactly what is missing and how to fix it
(verified: bad `OTILIB_ROOT` → exit 1, missing conda env → exit 1). All of
`OTILIB_ROOT`, `OTILIB_CONDA_ENV`, `CONDA_SH` are overridable.

**Correction to an earlier report.** An earlier acceptance summary stated
"OTILib is not available in Windows or WSL; the pass path could not be
exercised." **That was wrong.** The probe used WSL's *system* `python3` without
activating the `pyoti` conda env, so it missed a build that was present all
along (the extension is compiled for the env's Python 3.9, not system Python).
The WSL script exists precisely so this cannot happen again.

**On Windows (no WSL): the OTILib tests SKIP cleanly** — `otilib_available()` is
`False`, and `python -m pytest` reports 3 skips with a reason, never a silent
pass. OTILib does not build natively on Windows.

Two honest caveats:
- `pytest` had to be installed into the `pyoti` conda env (`python -m pip install
  pytest`); `scripts/setup_otilib.sh` does not install it. The script detects its
  absence and prints the exact install command rather than failing obscurely.
- The Python residual path (`resasm run` with `residual.type: python`) hard-requires
  OTILib in `oti_global.solve_python`, so it cannot run on Windows regardless of the
  `backend:` setting — the advertised `dual1` backend is currently unreachable there.
  Pre-existing; not changed.

---

## Part II — Crystal-plasticity C3D8+UMAT backend (original Step-10 criteria)

The following is unchanged; it is now the first verified backend
(`formulations/solid_c3d8_finite_strain.py` + `materials/crystal_plasticity_adapter.py`),
with tests under `tests/cp_c3d8_umat/`.

## Environment constraint (root cause of what is / isn't verifiable here)

**Abaqus is not installed in this environment** (`which abaqus` → none). Python
3.12 + numpy and gfortran are available. Consequently:

- Everything that does **not** need Abaqus is built AND verified here to machine
  precision (manifest, parsing, the finite-element residual math, the tangent
  machinery, the UMAT-driver plumbing).
- Everything that **needs Abaqus** (running the job, ODB export, comparing to
  reactions and to `STATEV`/`S` history, compiling the real UMAT with ifort) is
  built **ready-to-run and documented**, but is **not executed** — and is
  labelled ⏳, never ✅.

## How this was built (process)

A shared [`residual_core/CONTRACT.md`](residual_core/CONTRACT.md) fixed the
conventions first. Then **parallel agents** built the modules (manifest, parser,
residual/tangent, Abaqus IO + Fortran adapter); a second wave of **adversarial
cross-audit agents** each reviewed modules they did *not* write, re-deriving the
math from scratch and re-running the builds; confirmed findings were fixed and
re-verified. Every deliverable was also re-run by the orchestrator independently
(not trusting agent self-reports).

---

## Step-10 criteria

### 1. The manifest correctly classifies the downloaded examples — ✅ MET
- All **105** `.inp` files under `sources/{permissive,copyleft,license-unknown}`
  inventoried in `residual_core/example_manifest.{md,json}`.
- Groups: **A=4** (C3D8-only CP UMAT — the target), **B=3** (non-C3D8 continuum),
  **C=29** (UEL/cohesive/`U1`), **D=56** (all copyleft + license-unknown,
  reference-only), **fragment=13** (`*Include` mesh pieces).
- Independently re-verified by audit: counts match a fresh recount and on-disk
  `find`; the **license invariant holds exactly** (all 56 copyleft/unknown files
  are D — none can leak into the permissive core); the 4 priority files are all
  group A with correct element type / user-material / constants / SDV.

### 2. At least one permissive C3D8 + UMAT CP example is parsed — ✅ MET
- `abaqus_inp_parser.py` parses `HCPnoTwin/Compression111.inp`: **216 nodes, 125
  C3D8 elements, material CPuranium (user-material, 11 constants, 125 SDV)**, all
  BCs, resolving the tricky part/assembly `Set-1` name collision correctly.
  Independently re-derived from the raw file by audit — every number matches.
- Also parses the larger group-A jobs (DiscreteTwin 8035 elems, PyCiGen 2940
  elems/95 grains) without error.

### 3. Abaqus ODB fields are exported — ⏳ BUILT, AWAITING ABAQUS
- `extract_abaqus_fields.py` (Abaqus/Python-2.7) emits the CONTRACT-§5
  `fields.json` (undeformed nodes, connectivity, per-frame `U`/`RF`, IP `S` in
  Abaqus order with IPs sorted, `SDV`). odbAccess API usage audited and
  **exercised against a hand-built fake ODB** (grouping/sorting/schema correct);
  `py_compile` clean; 2.7-safe. Cannot run for real without Abaqus.

### 4. Stress-driven residual assembly reproduces Abaqus RF / free residuals — ⏳ PARTIAL
- **The assembly math is verified offline to machine precision** (this is the
  substance of the check, minus Abaqus' actual numbers):
  - Divergence-theorem patch test, real Compression111 element + a distorted
    hex, random + uniaxial stress: `∫BᵀσdV` = surface-traction consistent forces
    to **5e-16**; self-equilibrium to 1e-13.
  - **Linear-stress** patch + global body-load test (catches per-IP weights and
    intra-element redistribution that a uniform field can't): **1e-13**.
  - Finite-strain force is **frame-objective to 5e-16** and reduces to
    small-strain at `u=0` exactly — i.e. it really integrates Cauchy stress over
    the current config, the Abaqus nlgeom convention.
- `stress_driven_residual.py --fields fields.json` computes free-DOF residual and
  compares prescribed-DOF internal force to `RF` (±sign) — **ready**, but the
  actual RF comparison ⏳ awaits an ODB.
- Caveat (documented): the Abaqus **IP ordering** is not offline-verifiable and
  must be confirmed against the first ODB before trusting a non-uniform state.

### 5. UMAT replay reproduces Abaqus STRESS and STATEV history — ⏳ PLUMBING PROVEN, REAL UMAT PENDING ifort
- `umat_adapter_fortran/` marches a full DFGRD history in one process, persisting
  `STATEV` + the `/UMPS/` common block across increments (never a final-step
  jump). The **37-arg UMAT signature matches `umat.for` exactly** and the STATEV
  layout claims match `kmat.f` (both line-by-line audited).
- Verified with a **mock** elastic UMAT: 5-increment replay, `STATEV(35)`
  accumulates `0→1e-5→…→1e-4` (proves history dependence), identity increment →
  zero stress. This proves the driver/orchestrator plumbing.
- The **real** Grilli UMAT does not compile under gfortran (Cray-pointer/`target`
  twin arrays + an ifort `trace()` kind mismatch — confirmed, exactly two
  errors, nothing hidden). It requires **Intel ifort + Abaqus (MKL)**; the build
  route is documented. So STRESS/STATEV-vs-ODB ⏳ awaits that toolchain + a run.

### 6. Residual from replayed UMAT stress matches Abaqus RF / free residuals — ⏳ AWAITING 4 & 5
- This is the composition of (5)→(4): once the real UMAT replay produces per-IP
  stress and the ODB exists, feed that stress through the (already verified)
  assembler and compare to RF. Both halves are built; the end-to-end number
  awaits Abaqus + ifort.

### 7. Tangent check passes or clearly reports the DDSDDE mapping error — ✅ MET
- Small-strain element tangent `K = ∫BᵀDB` verified by finite differences:
  relative (Frobenius) error **1.6e-16**.
- Finite-strain force Jacobian at fixed σ verified against an exact analytic
  Jacobian (**5e-11**), and the audit confirmed this object is distinct from the
  conventional geometric stiffness (a docstring that conflated them was fixed).
- The report **clearly states the still-open mapping**: turning the *finite-strain
  UMAT `DDSDDE`* into the finite-strain element material tangent `𝔻` (objective
  rate + geometric split) to match Abaqus `AMATRX` is deferred to the Abaqus
  comparison — see `tests/c3d8_tangent_fd_check/README.md`. Numbers in
  `tests/c3d8_tangent_fd_check/last_report.txt`.

---

## Summary

| # | Criterion | State |
|---|-----------|-------|
| 1 | Manifest classification | ✅ verified |
| 2 | Parse a C3D8 CP UMAT example | ✅ verified |
| 3 | Export Abaqus ODB fields | ⏳ built, needs Abaqus |
| 4 | Stress-driven residual = Abaqus RF | ⏳ math verified offline; RF compare needs ODB |
| 5 | UMAT replay = ODB STRESS/STATEV | ⏳ plumbing proven; real UMAT needs ifort+Abaqus |
| 6 | Residual from replayed stress = RF | ⏳ both halves built; needs 4+5 |
| 7 | Tangent check / mapping error reported | ✅ passing + open item stated |

**Fully met and verified: 1, 2, 7.** **Built, ready, and verified as far as is
possible without Abaqus: 3, 4, 5, 6** — the remaining gap for each is a real
Abaqus run (and, for 5/6, an ifort build of the MIT UMAT), not missing code.

## Explicitly NOT done (per the brief's guardrails)
No generic elastic MVP; no OTI/HYPAD; no UEL/cohesive path (inventoried
separately as group C); no new/expanded CP physics; no 150-parameter model; no
final-step parameter overloading — the replay marches the full increment history
because `STATEV` is history-dependent.

## Immediate next actions when an Abaqus license is available
1. Run `Compression111` with the UMAT; export `fields.json` (`--frames all`).
2. Confirm the C3D8 **IP ordering** with a single-element spatially-varying-stress job.
3. Run Mode-1 (`stress_driven_residual.py --fields …`) at converged frames → expect
   free-DOF residual ≈ 0, prescribed-DOF internal force = ±RF.
4. Build the real UMAT (ifort+Abaqus), run Mode-2 replay increment-by-increment,
   compare STRESS↔S and STATEV↔SDV, then close Mode-1 with replayed stress.
