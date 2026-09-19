# Review guide

A guided tour of Residual_Assembler for a reviewer: what the program is, where
each capability lives in the code, which tests and scripts verify it, and how
to rerun them. Read it top to bottom. Sections 1 to 4 cover the sensitivity
workflow users run; sections 5 to 11 cover the assembly framework underneath.
Most steps run offline on Python and NumPy; the ones that need Abaqus, a
Fortran compiler or OTILib say so.

For installation see [docs/INSTALL.md](docs/INSTALL.md); for what is supported
and how it was measured, [STATUS.md](STATUS.md).

## 1. What the program does

Residual_Assembler computes parameter sensitivities of a converged
finite-element analysis without re-running it. From `Analysis.inp`,
`Analysis.odb`, a compiled OTI material provider (`OTI_UMAT.obj` and
`Mapping.json`, built by the companion UMAT-OTI) and a
`sensitivity_request.json`, it replays the material at every integration point
and increment, assembles `R`, `K` and `dR/dp`, solves `K du/dp = -dR/dp` and
writes the requested sensitivities.

```bash
resasm request --model Analysis.inp --odb Analysis.odb \
    --material OTI_UMAT.obj --request sensitivity_request.json --out results
```

It also assembles the global residual from its ingredients,

```
R(y, q, a, t) = F_internal(y, q, a, t) - F_external(t) + F_constraints(y, t)
```

outside the solver, and differentiates residuals users supply directly or
through a black-box executable.

## 2. Where the sensitivity workflow lives

| Piece | Code |
|---|---|
| Command routing (bounded engine or history engine) | `residual_core/ui/cmd_request.py`, `residual_core/ui/cmd_history.py::route_request` |
| Bounded J2 engine | `residual_core/replay/presentation.py`, `presentation_inputs.py`, `request_reductions.py` |
| History replay engine | `residual_core/replay/history.py`, `history_material.py`, `history_inputs.py`, `history_outputs.py`, `history_verify.py` |
| ODB export (Abaqus Python) | `residual_core/io/abaqus_odb_export.py`, `residual_core/replay/odb_export_npz.py` |
| Provider link and ABI checks | `residual_core/replay/path_material.py`, `abi.py`, `objlink.py` |
| GUI (thin front end over the CLI) | `residual_core/app/`, launched by `scripts/app.py` |

The mathematics, tolerances and supported deck subset are in
[docs/REPLAY_HISTORY.md](docs/REPLAY_HISTORY.md); the request format and the
bounded engine in [docs/REQUEST_INTERFACE.md](docs/REQUEST_INTERFACE.md).

## 3. How the sensitivities are checked

Each check compares against a reference that does not share the code path
under test:

- **Against the ODB.** Replayed stress, state and reactions are compared with
  the ODB at every integration point and increment; any excess is a failure,
  never a warning.
- **Against finite differences of the ORIGINAL UMAT.** The untransformed
  routine compiled into the same provider object is re-equilibrated over the
  whole model at perturbed parameters (`resasm history --verify fd`,
  `resasm request --validate`).
- **Against Abaqus reruns** at perturbed parameters
  (`examples/cantilevers/run_fd.sh`, `fd_reference.py`).
- **Against Euler's homogeneity identity**, which needs no finite differences:
  for models homogeneous of degree one in their stress-dimensioned parameters,
  the parameter-weighted sensitivities of each stress and reaction sum to the
  value itself at every increment.
- **Against closed forms** for the bounded J2 request (analytic uniaxial
  derivatives) and the direct and black-box residual templates.

Tests: `tests/replay_history/` (history engine, provider, outputs, the
committed Abaqus beam), `tests/integration/test_presentation_request.py`
(bounded engine), `tests/integration/test_connected_j2.py` (connected replay),
`tests/gui/` (the Sensitivity Request screen, including a headless-browser run
selected with `-m gui`).

Every quantitative claim, with the command that reproduces it and the value
measured, is in [docs/VERIFICATION_RECORD.md](docs/VERIFICATION_RECORD.md);
`verification/run_all.py` reruns them.

## 4. Clean installation

`scripts/clean_install_gate.py` builds wheels of both repositories from clean
trees, installs them in a new virtual environment with no inherited Python
path, and runs the workflows from the installed commands only: provider build,
`resasm request` on a genuine ODB with an analytic check, the same request with
every read of a Fortran source denied, both GUIs up to HTTP readiness, and
(with `--cantilever`) the full-size J2 cantilever. The result from clean clones
of the published branches is in
[docs/evidence/final_clean_clone.md](docs/evidence/final_clean_clone.md).

## 5. The assembly framework

The core assembler contains no element or material mathematics: it dispatches
each element to a registered **formulation backend** and each material point
to a registered **material backend**, then scatters. See
[residual_core/docs/architecture.md](residual_core/docs/architecture.md).

The argument that it is model-agnostic, made from the code, is in
[residual_core/docs/minimality_proof.md](residual_core/docs/minimality_proof.md):
no physics in the core, heterogeneous per-node DOFs, backends that declare
their capabilities, and a requirements engine that asks for the minimum
missing input. The same assembler drives stress-driven solids, elastic replay,
crystal plasticity, truss, beam and a mixed model unchanged.

## 6. Implemented backends

Nine formulation backends and four materials (`resasm backends`):

- `stress_driven_c3d8` (Mode 1), `solid_c3d8_small_strain` and
  `solid_c3d8_finite_strain` (Mode 2), `uel_direct` (Mode 3, skeleton),
  `truss2`, `beam2`, `nonlinear_spring1` and `nonlinear_bar1` (Mode 4),
  `shell_placeholder` (contract only);
- materials: `isotropic_elastic`, `compressible_neo_hookean`, `umat`,
  `crystal_plasticity`.

Crystal plasticity is one backend, not the purpose of the framework. Full
status: [residual_core/docs/support_matrix.md](residual_core/docs/support_matrix.md).

## 7. External UMAT and UEL cases

Public Abaqus user subroutines were catalogued across diverse physics
(viscoelasticity, J2 plasticity, damage, hyperelasticity, cohesive laws,
coupled chemo-mechanics and crystal plasticity). No external code is vendored;
copyleft, licence-unknown and other-solver sources are reference-only.

- Provenance and licences: [sources/external_manifest.md](sources/external_manifest.md),
  [sources/external_manifest.json](sources/external_manifest.json),
  [sources/LICENSE_AUDIT.md](sources/LICENSE_AUDIT.md).
- Verification cards: [tests/verification_zoo/](tests/verification_zoo/), one
  `CASE.md` per case with source, licence, backend, required inputs, expected
  outputs, offline and Abaqus tests, status and the next missing item.

Each case has a **residual class** (A to H) and a **minimum verification mode**
(1 stress-driven, 2 material replay, 3 direct UEL, 4 native formulation); see
[tests/verification_zoo/README.md](tests/verification_zoo/README.md) and
[manifest_summary.md](tests/verification_zoo/manifest_summary.md).

## 8. What is verified offline

```bash
# one-shot guided demo
bash scripts/demo_framework.sh

# or the individual suites
python tests/framework/test_assembler.py               # core == kernel (bitwise)
python tests/framework/test_truss2_backend.py          # EA/L
python tests/framework/test_beam2_backend.py           # PL^3/3EI
python tests/framework/test_mixed_model_dispatch.py    # physics-blind dispatch
python tests/framework/test_dof_manager_mixed.py       # heterogeneous DOFs
python tests/framework/test_requirements_negative.py   # minimum missing input
python tests/framework/test_public_api.py              # facade only
python tests/framework/test_neutral_io.py              # neutral JSON round-trip
python -m residual_core.ui.examples                    # runnable API examples
```

A transcript of real output is in
[residual_core/docs/demo_transcript.md](residual_core/docs/demo_transcript.md).
The whole offline suite is
`python -m pytest -q -m "not abaqus and not arc and not network"`.

## 9. What uses Abaqus

Abaqus (tested: 2021.HF5) is needed only to run jobs and export ODBs. The
sensitivity computation never calls it, and tests marked `abaqus` are the only
ones that need it. The replay engines compare against the ODB whenever one is
supplied. The earlier ready-to-run comparison scripts for the assembly
framework remain available and skip cleanly without Abaqus:

```bash
python scripts/run_abaqus_validation.py --job J --inp model.inp --user umat.for
abaqus python scripts/extract_odb_fields.py --odb J.odb --out fields.json
python scripts/compare_residuals.py   --model model.json --fields fields.json
python scripts/compare_umat_replay.py --reference ref.json
python scripts/compare_uel_rhs.py     --uel-dump uel.json
```

The order in which the framework's remaining Abaqus comparisons close is in
[residual_core/docs/abaqus_validation_roadmap.md](residual_core/docs/abaqus_validation_roadmap.md);
plans: [umat_validation_plan.md](residual_core/docs/umat_validation_plan.md),
[uel_validation_plan.md](residual_core/docs/uel_validation_plan.md).

## 10. How to run demos and tests

- **Demo:** `bash scripts/demo_framework.sh` (offline).
- **API:** `from residual_core import ResidualProblem`; see
  [residual_core/docs/user_interface.md](residual_core/docs/user_interface.md).
- **CLI:** `resasm --help` (or `python -m residual_core.ui.cli --help` without
  an installed console script); every command is in
  [docs/CLI_GUIDE.md](docs/CLI_GUIDE.md).
- **Tests:** `python -m pytest -q -m "not abaqus and not arc and not network"`;
  the `tests/framework/*.py` scripts above also run on their own and print
  `PASS` or `FAIL`.
- **Repository standards:** `python tools/audit_repository_standards.py`.

## 11. Guardrails

OTILib is the hypercomplex backend: an external dependency, reported cleanly
when absent, with no silent fallback. Dual1 is only a legacy first-order smoke
test. No imitation OTI library is bundled, and no unsupported physics or
restricted or licence-unknown code is integrated. Crystal plasticity remains
one verified example backend among many.
