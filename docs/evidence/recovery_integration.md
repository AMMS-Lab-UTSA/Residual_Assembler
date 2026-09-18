# Recovery Integration Evidence

Date: 2026-09-18. Scope: W2 plus bounded W5, only `imq-ra-recovery` implementation
changes. Read `../imq_BRIEF.md`, branch audit, W1 evidence, W3 provider docs and
W3 recovery evidence before implementation. No subagents, commits, pushes,
branch switches, detached commands, process kills, privilege prompts, corpus
runs or Abaqus jobs. Old Claude worktrees were read only. Existing W1 kernels,
`resasm_user/__init__.py`, reporting edits and existing tests were preserved.
Only one CLI registration line was added to the already modified CLI file.

## Recovery Provenance

Inspected `../imq-ra-W2-replay/residual_core/replay/path_material.py`,
`core/field_sensitivity.py`, `replay/kinematics.py`, runtime exports and tests.
The first two compared byte-identically with their source-history counterparts.
Restored absent replay/runtime/interface files, derivative-field loader, field
tests and reference example from `origin/cross-platform-hardening`, commit
`bfde4d027652b1b60d306a42bdc6746f64ddb24f`, with `git archive | tar --keep-old-files`.
No existing recovery file was overwritten. Recovered the needed B-bar operators
from the inspected untracked W2 kinematics via apply_patch, retaining them under
replay ownership rather than modifying the W1 kernel. New behavior/tests/docs
were all edited with apply_patch.

Real Abaqus evidence restored from the same committed upstream history:
`tests/abaqus_derivative_export/nonuniform_c3d8.inp`, verification elastic UMAT,
SDV layout, derivative fields, IP stresses, nominal U and four perturbed U files.
Their provenance is `nonuniform_c3d8.odb`, Abaqus 2021, Step-1/frame 1. These are
real elastic exports, not generated J2 data. Only this elastic subset is used;
the retained historical README also describes upstream files outside this slice.

## Implementation and Checks

- Fresh W3 transformation/compilation; PathMaterial EVAL/MARCH tested without
  mocks or old prebuilt objects. Dimensions `(6,4,1,4)`, parameter order
  `E,nu,SIGY0,H`, 18 EVAL arguments, 11 MARCH arguments, physical state separate
  from derivative buffers. Full object SHA-256 checked when declared.
- Whole-material ORIGINAL FD verifies stress/state at three steps through
  elastic/plastic loading, unloading and reverse plasticity; EVAL/MARCH agree.
- Fingerprint-gated strain-history seeding through incoming stress derivatives
  carries both stress and state sensitivities, verified against changing-strain
  ORIGINAL paths. No J2 material-law rewrite or constitutive FD in the solver.
- Nonlinear multi-IP C3D8 assembly and constrained Newton with line search.
  Every trial starts from committed history. After convergence the causal
  sensitivity solve carries previous equilibrium derivatives into the next step.
- Independent whole-model ORIGINAL FD verifies all four parameters for U,
  stress and state; separate fixed-path material FD, all-IP branch checks and
  last-two-step plateau checks. Fixed-path local solve is explicitly not du/dp.
- CLI solve, record reload/replay, negative input and failure summary tests.
  GUI AppTest runs the real compiled CLI and checks the verified output.
- Regression checks reject unsupported physics, mismatched source provenance,
  corrupted recorded physical state, invalid path lengths/time increments, and
  live OTILib scalars at extracted-real field/ABI boundaries.
- Reproducer tests fresh builds, both commit hashes, full object hash, plot/CSV
  output and nonzero exit with a failed manifest for missing input.

The first restoration run had 17 passes/2 failures (missing local import/assets),
then 18/1 (remaining schema dependency), then 19/0. These were resolved by
restoring actual dependencies, not changing test expectations. The initial
whole-model check exposed Newton unloading overshoot and precision noise:
added Armijo backtracking and tightened primal convergence from `1e-11` to
`1e-14`. The derivative acceptance limit stayed `2e-6`. One GUI test harness
failure was fixed by using the existing state initializer and `btn_` key.

## Exact Environment and Commands

All commands synchronous, from `/home/ammslab3/softwarex_work/imq-ra-recovery`.
Use the following environment for the commands below (actual runs used these
same assignments inline with `env`):

```sh
cd /home/ammslab3/softwarex_work/imq-ra-recovery
export PY=/home/ammslab3/softwarex_work/.venv/bin/python
export PYTHONPATH=/home/ammslab3/softwarex_work/imq-ra-recovery:/home/ammslab3/softwarex_work/imq-umat-recovery/src:/home/ammslab3/otilib/build_py311
export UMAT_OTI_REPO=/home/ammslab3/softwarex_work/imq-umat-recovery
export PYOTI_PATH=/home/ammslab3/otilib/build_py311
export OTILIB_ROOT=/home/ammslab3/otilib/build_py311
export RUN_OTILIB_TESTS=1

env GIT_TERMINAL_PROMPT=0 bash scripts/init_permissive_sources.sh --required-only

"$PY" -m pytest -q tests/integration/test_connected_j2.py \
  tests/framework/test_field_sensitivity.py tests/framework/test_binary_compat.py \
  tests/framework/test_interface_versions.py tests/framework/test_user_layer.py \
  --junitxml=.pytest_cache/recovery_integration_focused_final.xml

"$PY" scripts/reproduce_imqcam_pipeline.py --skip-abaqus \
  --out .pytest_cache/recovery_integration_final

"$PY" -c 'import os,sys,pytest,residual_core,umat_oti,pyoti; assert os.getcwd()=="/home/ammslab3/softwarex_work/imq-ra-recovery"; assert residual_core.__file__.startswith(os.getcwd()); assert umat_oti.__file__.startswith(os.environ["UMAT_OTI_REPO"]); print("RECOVERY IMPORTS",sys.executable,residual_core.__file__,umat_oti.__file__,pyoti.__file__,flush=True); sys.exit(pytest.main(["-q","-ra","-m","not abaqus and not arc and not network","--junitxml=.pytest_cache/recovery_integration_offline_final.xml"]))'
```

For another reproduction use a fresh output name; existing output directories
are deliberately refused. No generated output listed here is a required input.
The provider and both public CLI subprocess invocations, including absolute
paths and exit codes, are recorded in `private/manifest.json` and command logs.
The equivalent standalone commands are in `../IMQCAM_WORKFLOW.md`.

Verified imports: recovery `residual_core`, recovery `umat_oti/src`, and
`/home/ammslab3/otilib/build_py311/pyoti`. Python 3.11.7; NumPy 2.4.6;
matplotlib 3.11.1; pytest 9.1.1; Streamlit 1.62.0; GNU Fortran 9.4.0.
`abaqus` availability: `/usr/bin/abaqus`; no licensed job was needed or run.
Initialized only the documented MIT/BSD-tier Oxford dependency at
`85102ed35dc4592edc5fd4aaec542fa63e7162fd`. Restricted source tiers remained empty.

## Results

| Gate | Passed | Failed | Skipped |
| --- | ---: | ---: | ---: |
| Restored field/runtime/interface tests | 19 | 0 | 0 |
| Final focused (including fresh reproduction and GUI) | 47 | 0 | 0 |
| First full offline, before dependency/environment repair | 334 | 0 | 19 |
| Final full offline (354 cases, 41.737 s) | 343 | 0 | 11 |
| Genuine OTILib subset of final full suite | 15 | 0 | 0 |

JUnit was parsed after completion: zero failures/errors, exactly 15 existing
`test_otilib_*` cases, none skipped. No tests were disabled or marked skip/xfail.
All new integration tests ran. Editor diagnostics found no errors in the new
integration modules, field solver, CLI/GUI, script and tests.

Final retained reproduction: `.pytest_cache/recovery_integration_final`.
`manifest.json` has `passed=true`. Measured worst scaled errors over all four
parameters and all FD step sizes:

| Check | Maximum |
| --- | ---: |
| Free equilibrium residual / max(external force, 1) | 2.842170943040401e-15 |
| ORIGINAL versus lifted nominal U | 1.5978127228782124e-15 |
| Total-history U derivative versus ORIGINAL FE FD | 1.0247424053258358e-8 |
| Total-history stress derivative versus ORIGINAL FE FD | 9.720407922183769e-7 |
| Total-history state derivative versus ORIGINAL FE FD | 2.8375776287807234e-8 |
| Fixed-path stress derivative versus ORIGINAL path FD | 2.184842165358705e-7 |
| Fixed-path state derivative versus ORIGINAL path FD | 9.95439789547374e-9 |
| Worst FD plateau discrepancy (total stress) | 1.0284862393444268e-6 |

Every derivative and plateau error is below the unchanged `2e-6` limit.
The maximum absolute difference between true du/dp and the fixed-path local
solve is `0.0012843052410761524`, proving why their semantics must stay separate.

Real archived elastic Abaqus checks: du/dE error `2.0823535819322397e-5`
(limit `3e-5`), du/dnu `5.454066796529612e-6` (limit `1e-5`), IP stress
`7.647495682963444e-8` (limit `1e-4`), equilibrium `1.5715691814577325e-8`
(limit `1e-5`). Deliberately permuted IPs produce equilibrium error
`0.4752019982761178`, correctly exceeding the `0.1` rejection threshold.
These separate fixture tolerances reflect the original float32 ODB FD floor;
they are not the stricter double-precision J2 tolerances.

## Fingerprints

Git HEAD alone does not describe the dirty recovery state. The manifest records
per-source SHA-256 hashes and aggregate hashes in addition to both HEADs:

```text
RA HEAD:   69ac77cf4865ef5cfc2a46fad987bb7fb6d73e21
UMAT HEAD: d34c7fff51e13285efaf13a57917fa34034f57cb
RA source tree:   cec8f09292eab723312703e9cdf8803e7570daa442ddaa2c288ab9d43a6f512d
UMAT source tree: 90bdb927e50d83b9e2b1afdca63a8dc292555454f11250c6fd5e2a74a17a7f6f
ORIGINAL source: 9b779f0c6cadf9c493068d84cd18bbf91f91ae14cb08aa0e7c8ac4ca119a9b40
fresh object:    87c38664872e1d90d825c301b07ab2f46b95df28e5b8906432c664a1e8f11403
fresh contract:  b59093a4250f7dbe3dbb7bc7061f4b185eedec482d7e90fd2b1624f01628c402
```

Hashes describe the retained run, not a promise of binary identity across
compiler/build-directory environments. The build sources and command logs are
retained for diagnosis. This manifest predates the documentation-only additions,
including the example README; its per-file list identifies exactly what was
hashed at execution time. The executed Python/model inputs did not change after
the retained run.

## Remaining Boundaries

Eleven existing cross-reader checks cannot run because the actual frozen
`corpus_run/pass11/results/store_verification.jsonl` is absent. They are NOT
passes; no data was fabricated, substituted from another corpus pass, or fetched
through a forbidden whole-corpus run. The first run's producer-path skips were
fixed with `UMAT_OTI_REPO`, and all eight missing-source skips were resolved.

No fresh J2 Abaqus comparison, full-sized J2/FCC cantilever, generic path UMAT,
finite-strain physics, arbitrary BC/load/geometry derivatives, higher order,
performance scaling or clean-install verification is claimed. The single-cube
cyclic loads deliberately amplify history and can reach about 10% displacement;
this validates small-strain equations, not large-deformation physical accuracy.
Legacy C-package/job modules are recovered dependencies, not newly validated
generic J2 interfaces; the supported connected path is `resasm replay` through
PathMaterial. No completion claim is made for the full 274-row program.