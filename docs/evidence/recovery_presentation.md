# Recovery Presentation Evidence

Date: 2026-09-18. Implementation only in `imq-ra-recovery`; test artifacts under
`../imq_abaqus/recovery_presentation/imqrp_reference`. All pre-existing dirty
changes preserved. No commits, staging, pushes, branches, nested agents or
original development worktree edits. One invocation controlled terminal execution;
all commands synchronous. Browser harnesses owned/terminated only their exact
server PIDs. No server is left running. No source was deleted.

## Grounding And Changes

Read shared brief, branch audit, workflow and recovery integration evidence.
Read slide text and original PPTX XML for slides 9-12 and 14-18. PPTX SHA-256:
`55688498a351ea403788b49faf7a7ac26758c84483588fa54dd3e95f0e8ae81d`.
Slide 12 names are the CLI/primary GUI inputs and public output names.
Slide 11's generated mapping is discovered or explicitly selected, never
manually reconstructed. FCC/full-size presentation capabilities are not claimed.

New: `replay/presentation.py`, `presentation_inputs.py`, `request_reductions.py`,
`ui/cmd_request.py`, example INP/request, presentation integration tests,
genuine exported-field fixture, proof/browser scripts, interface guide.
Modified: one CLI registration, primary GUI request screen, existing exporter
strict mode, connected replay's explicitly optional ODB precision checks.
The existing exporter was already present; no history recovery was needed.

The GUI obeys the existing thin-CLI architecture: GUI `_run` -> `cli.main` ->
`cmd_request` -> shared `presentation.run_request`. No alternate mathematics.
Plastic history never uses legacy `request_contract.replay_sensitivities`.
Engine arrays are checked as `du_dp(ndof,4)`, `dsigma_dp(nelem,8,6,4)` and
`dstate_dp(nelem,8,1,4)`. RF derivatives use `Rp_history + K @ du_dp`.
LAST/list output filtering occurs only after all earlier history was replayed.

Mapping full object hash, source fingerprint, dimensions/ABI, PROPS/direction
order and derivative/Voigt layouts are validated. Completed mapping producer:
`../imq-umat-recovery/src/umat_oti/provider/build.py`. No generator edits.
Alias `Mapping.json` is the unchanged generated contract; conflicting sidecars
fail unless explicitly selected. Only compiled provider and mapping are needed.

## Exact Environment And Commands

Run from `/home/ammslab3/softwarex_work/imq-ra-recovery`:

```sh
export PY=/home/ammslab3/softwarex_work/.venv/bin/python
export PYTHONPATH=/home/ammslab3/softwarex_work/imq-ra-recovery:/home/ammslab3/softwarex_work/imq-umat-recovery/src:/home/ammslab3/otilib/build_py311
export UMAT_OTI_REPO=/home/ammslab3/softwarex_work/imq-umat-recovery
export PYOTI_PATH=/home/ammslab3/otilib/build_py311
export OTILIB_ROOT=/home/ammslab3/otilib/build_py311
export RUN_OTILIB_TESTS=1

"$PY" -m pytest -q tests/integration/test_presentation_request.py \
  --junitxml=.pytest_cache/presentation_final_focused.xml
"$PY" -m pytest -q tests/integration/test_presentation_request.py \
  tests/framework/test_gui_is_a_thin_cli_front_end.py \
  --junitxml=.pytest_cache/presentation_gui_architecture.xml
"$PY" -m pytest -q -ra -m 'not abaqus and not arc and not network' \
  --junitxml=.pytest_cache/presentation_offline_final.xml
```

Verified Python 3.11.7 and imports from recovery `residual_core`, recovery
`umat_oti/src`, and `/home/ammslab3/otilib/build_py311/pyoti`. Provider compiled
fresh using gfortran; no core math mocks. Test-only exporter substitution is
explicitly named `offline-test-transport`, never counted as actual ODB evidence.

| Gate | Passed | Failed | Skipped | Duration |
| --- | ---: | ---: | ---: | ---: |
| First scalar edit, immediately checked | 8 | 0 | 0 | 0.17 s |
| Connected history plus initial reductions | 28 | 0 | 0 | 25.55 s |
| Final presentation regressions | 41 | 0 | 0 | 8.54 s |
| Presentation plus existing GUI architecture guards | 49 | 0 | 0 | 8.47 s |
| Final full RA offline suite | 384 | 0 | 11 | 51.49 s |
| Genuine OTILib subset of full suite | 15 | 0 | 0 | included |

JUnit parsed explicitly: 395 cases, 11 skipped, zero failures/errors. Eleven
existing skips require absent `corpus_run/pass11/results/store_verification.jsonl`.
They are not passes; no replacement corpus was fabricated. Existing connected
J2 ORIGINAL whole-model/history FD and genuine OTILib regressions ran in full.

The first full run reported 383 passes/1 failure/11 skips: direct service import
in the GUI violated its architecture guard. GUI was routed through its existing
CLI bridge; the guard was not changed. One subsequent GUI layout edit had a
syntax error caught immediately by its focused test, then fixed and retested.
All final executable files were covered by the passing full run.

## One Genuine Licensed Job

Inspected `/usr/bin/abaqus` and
`/opt/intel/oneapi/compiler/latest/linux/bin/intel64/ifort` before launch.
Only one analysis job was launched, prefix `imqrp_`:

```sh
source /opt/intel/oneapi/compiler/latest/env/vars.sh
"$PY" scripts/reproduce_presentation_request.py prepare \
  --work /home/ammslab3/softwarex_work/imq_abaqus/recovery_presentation/imqrp_reference
```

The script synchronously executed:

```text
/usr/bin/abaqus job=imqrp_j2 input=/home/ammslab3/softwarex_work/imq-ra-recovery/examples/presentation_request/Analysis.inp user=/home/ammslab3/softwarex_work/imq-umat-recovery/parameter_sensitivity/models/m3_j2/umat.for cpus=1 interactive
```

Working directory:
`../imq_abaqus/recovery_presentation/imqrp_reference/developer/imqrp_j2`.
Abaqus 2021.HF5, ifort 2021.10.0; `Abaqus JOB imqrp_j2 COMPLETED`, exit 0;
`.sta` confirms `THE ANALYSIS HAS COMPLETED SUCCESSFULLY`.
Four converged increments, eight nodes, one C3D8, eight IPs, final plastic
equivalent strain about 0.025. This is a developer/test-only reference job,
not a collaborator requirement. No second job was used to fix import/replay.

Retained `developer/abaqus.log`, `developer/manifest.json`, job files and fresh
provider build. The manifest records the command, compiler, original source
hash and provider paths. The collaborator directory contains only:

```text
Analysis.inp
Analysis.odb
OTI_UMAT.obj
Mapping.json
sensitivity_request.json
```

## Source-Denied Collaborator Proof

```sh
"$PY" scripts/reproduce_presentation_request.py consume \
  --work /home/ammslab3/softwarex_work/imq_abaqus/recovery_presentation/imqrp_reference \
  --out results_final
```

The child runs from the isolated collaborator directory with explicit recovery
module paths. An audit hook rejects private `.for` reads and material-model
source directories. Only generated `private/link/path_shim.for` is exempt.
It is an audited execution check, not an OS security sandbox. The guard's
intentional blocked-open probe prints `PRIVATE SOURCE READ DENIAL ACTIVE`.
Neither source nor transformation is supplied to the request. The bundled
ORIGINAL binary symbol supplies the virgin tangent; that is not source access.

Underlying normal CLI command:

```sh
"$PY" -m residual_core.ui.cli request --model Analysis.inp --odb Analysis.odb \
  --material OTI_UMAT.obj --request sensitivity_request.json --out results_final
```

Actual result: `request executed: 4 scalar results; verified=False`, exit 0.
Private export command uses `/usr/bin/abaqus python` with `--frames all --strict yes`.
No solver job, synthetic seed, FD or production rerun occurred in consume.
All four increments replayed before LAST output selection. Three exact public
files and private fields/record/result were written. Public scalar results
include displacement, reaction, stress and state for four parameters.

Independent developer analytic check: homogeneous uniaxial plastic J2,
`U1=300/E+(300-SIGY0)/H`; zero-reference derivatives use absolute error.

| Measured check | Error |
| --- | ---: |
| Maximum scaled free residual | 4.120058472760017e-7 (limit 1e-5) |
| dU1/dE relative error | 1.6568627353483993e-7 |
| dU1/dSIGY0 relative error | 8.673617379884035e-16 |
| dU1/dH relative error | 8.839354487082118e-9 |
| dU1/dnu absolute error (analytic zero) | 1.1729503437354757e-9 (limit 1e-8) |
| Total support RF1 absolute error | 1.2360175412595709e-4 |
| S11 absolute error | 1.2360175412595709e-4 |
| EQPLAS absolute error | 2.209838492750471e-10 |

The ordinary public result correctly remains `verified=false`: execution
checks are not independent derivative verification. The analytic proof lives
separately in `analytic_checks.json`, `passed=true`, beside retained
`collaborator_invocation.log` and `collaborator_command.json`.

Initial proof attempts remain on disk: the source guard first blocked its own
generated shim; the exact shim was exempted. Actual ODB exposed fixed-DOF
residue about 1.9e-35 and float32 U reconstruction error in nominally zero
transverse stress about 1.23e-4. ODB-only BC and field-scaled tolerance handling
was added with regressions rejecting meaningful corruption. Synthetic solver
tolerances did not change. An initial analytic report divided a zero-reference
derivative by an arbitrary tiny scale; this was replaced with explicit absolute
error for zero references. These were failures/repairs, not claimed passes.

## Browser Evidence

```sh
"$PY" scripts/check_presentation_browser.py \
  --inputs /home/ammslab3/softwarex_work/imq_abaqus/recovery_presentation/imqrp_reference/collaborator \
  --out .pytest_cache/presentation_browser_cli_final
```

Real Playwright Chromium, actual ODB/provider, no transport mock. Selected
primary tab, filled four local file paths, executed the request through the
CLI, downloaded each artifact and checked exact filenames. Desktop 1440x1100
and mobile 390x844 screenshots; no document horizontal overflow or Streamlit
exception. The first mobile screenshot exposed the legacy sidebar obscuring
the form; default was changed to collapsed and screenshots rerun.

Final browser report: passed=true, server PID 1023122, temporary URL
`http://127.0.0.1:46425`, server_stopped=true. That URL is **not live**.
Harness terminates/waits that PID in finally, with exact-PID kill only as a
timeout fallback. Earlier browser PIDs 1009462 and 1010361 also stopped normally.
Retained screenshots, downloads, logs and JSON report are under the output above.

## Fingerprints And Remaining Limits

| Artifact | SHA-256 |
| --- | --- |
| Original UMAT | `9b779f0c6cadf9c493068d84cd18bbf91f91ae14cb08aa0e7c8ac4ca119a9b40` |
| Analysis.inp | `1cde25bd934403a4c3b341d7deca8d5a47897ae700f4329ab448d72d995b74d3` |
| Analysis.odb | `54745b97b515de159c3251e8a3ba0ed89ea1c65452e388ef34c8c81794e3b736` |
| OTI_UMAT.obj | `78b6632db73bc061fe241a0290e5aa6a47dac32089abd65e1483fa4812d3f71c` |
| Mapping.json | `acb21be0cfdc66c35acd3eb6735dd44e5da2077829e882231293c4c4aaf8c204` |
| Retained real fields.json | `98a8e3bc2e4157135033ee022961e71d288f13d86ce63460c75da103be85f950` |

Only 9,274-byte genuine fields and provenance are retained in the repository
fixture; no ODB/provider binary was added. Fresh builds need not be bit-identical.
Current scope is one small-strain static homogeneous pinned J2 C3D8/B-bar step,
virgin state, zero BCs and ramped nodal loads with every saved increment.
Dense scaling, float32 ODB error and documented unit-dependent tolerance floors
remain. Generic UMATs, finite strain, FCC/full-size cantilevers, pressure/contact,
amplitudes/multistep, nonzero BCs/initial state, higher derivatives, volume-weighted
reductions and von Mises output are not supported. Strict rejection is deliberate.
Explicit FD on float32 ODB may fail the older double-precision validation gate;
it is not enabled or claimed verified by normal use. No clean-install claim.