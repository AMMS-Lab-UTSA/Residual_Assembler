# Residual Assembler: Current Usage

Audited 2026-09-18 in the recovery working tree on Linux, Python 3.11.7,
gfortran 9.4.0 and genuine compiled OTILib. **Five examples ran and passed their
bounded numerical checks. This is not completion of the 274 requirements.**
No final-branch clean clone was tested. No Abaqus analysis was launched in this
audit; an existing genuine ODB was read using licensed Abaqus Python.

## What Works Now

Use the direct residual templates for scalar-generic order-two sensitivities;
the C3D8 stress-driven path for supplied integration-point stresses; the pinned
J2 replay for first-order total-history material sensitivities; and the bounded
neo-Hookean path for finite-strain assembly and first-order parameter sensitivity.
These are different supported paths, not a generic arbitrary-material FE engine.

The latest retained full offline suite is **419 passed, 11 existing skips, zero
failures/errors** ([fixture evidence](evidence/recovery_fixtures.md)). Historical
fixture arithmetic is not current material verification: only elasticity and
bundled J2 have current regenerated fixtures. Older finite-strain evidence's
suite failures were subsequently repaired without relabelling archived data.

## Installation

For the connected workflow use healthy Python >=3.10 with ctypes, ssl and venv.
RA alone declares Python >=3.9, but its UMAT companion requires >=3.10.
Linux/gfortran is the verified platform; Windows/macOS binary-provider support
is not established. RA requires NumPy; select `gui`, `yaml`, `test` extras for
these workflows. The companion's `paper` extra supplies plotting for the joint
reproducer. Do not select RA's historical `bridge` pin for this recovery pair.

With source checkouts available, the following is an installation recipe, not
a claim that a clean final-branch clone was tested here:

```sh
RA="$HOME/softwarex_work/imq-ra-recovery"
UMAT="$HOME/softwarex_work/imq-umat-recovery"
BASE_PYTHON="$HOME/anaconda3/bin/python3.11"
ENV=$(mktemp -d /tmp/resasm-env-XXXXXX)
"$BASE_PYTHON" -m venv "$ENV"
"$ENV/bin/python" -m pip install "$RA[gui,yaml,test]" "$UMAT[test,paper]"
"$ENV/bin/python" -m pip check
"$ENV/bin/resasm" --help
```

The previous [working-tree wheel gate](evidence/recovery_install.md) used a new
venv and scratch HOME, verified installed resources and the artifact-only J2
workflow. It predates later edits and is **not final-branch clean-clone proof**.
Examples/templates and some source discovery need the checkout, not just a wheel.
After setting the shared environment below, the joint wheel gate can be rerun
with a new external work directory:

```sh
"$PY" scripts/clean_install_gate.py --umat-repo "$UMAT" \
  --python "$BASE_PYTHON" \
  --odb "$RA/../imq_abaqus/recovery_presentation/imqrp_reference/collaborator/Analysis.odb"
```

The shared verification environment already exists. Use these exact selectors
for the remaining source-tree commands; the globally installed editable console
scripts can otherwise resolve to the original checkouts:

```sh
WORKSPACE="$HOME/softwarex_work"
RA="$WORKSPACE/imq-ra-recovery"
UMAT="$WORKSPACE/imq-umat-recovery"
PY="$WORKSPACE/.venv/bin/python"
BASE_PYTHON="$HOME/anaconda3/bin/python3.11"
export PATH="$WORKSPACE/.venv/bin:$PATH"
export PYTHONPATH="$RA:$UMAT/src:$HOME/otilib/build_py311"
export UMAT_OTI_REPO="$UMAT"
export PYOTI_PATH="$HOME/otilib/build_py311"
export OTILIB_ROOT="$HOME/otilib/build_py311"
export RUN_OTILIB_TESTS=1
cd "$RA"
resasm() { "$PY" -m residual_core.ui.cli "$@"; }
"$PY" -c 'import residual_core,umat_oti,pyoti.sparse; print(residual_core.__file__); print(umat_oti.__file__); print(pyoti.sparse.__file__)'
```

Both OTILib variables must name the healthy build, not its old source package.
Never install the unrelated PyPI `pyoti`. The existing external build is reused,
not redistributed or rebuilt here. Source-backed Oxford inspection examples
additionally need `bash scripts/init_permissive_sources.sh --required-only`;
the five examples below do not require that download. Missing compiler, OTILib,
source or ODB dependencies must be resolved or recorded as unavailable, not skipped
and counted as a successful scientific check.

## CLI Reference

[Captured help](evidence/usage_help.json) contains the actual stdout/stderr,
argv, cwd and exit code for the root and **all 17 public subcommands**, plus
example/gate scripts and Streamlit. Every listed help invocation exited zero.
Help proves argument availability, not successful physics for arbitrary inputs.
`--config FILE` is a global option and goes before the subcommand.

| Command syntax after `resasm` | Meaning / important limit |
| --- | --- |
| `backends` | Registry declarations and limits, not universal readiness |
| `modes` | Lists assembly modes |
| `inspect MODEL --detail` | Mesh and backend selection |
| `inspect-model MODEL --solution U.npy --material FILE --param NAME` | User-input inspection; optional arguments depend on the mode |
| `requirements MODEL --mode MODE --fields FIELDS` | Minimum missing inputs |
| `assemble MODEL --mode stress-driven --fields FIELDS --out R.npy` | Assembles residual; `--tangent` only when available |
| `verify MODEL --fields FIELDS` | Residual/reaction checks, not constitutive derivative certification |
| `sensitivity MODEL --params PARAMS --out PREFIX` | Parameter package; also `--param`, `--mode`, `--order`, `--backend`; order is backend-limited |
| `init --template python --out DIR` | Templates: python, blackbox, blackbox-order2, cpp, fortran; omit flags for wizard |
| `init-assembly --model MODEL --solution U.npy --material FILE --param NAME --out CONFIG` | Creates an assembly configuration, not a solve |
| `check CONFIG` | Configuration readiness |
| `run CONFIG` | Configured sensitivity workflow |
| `report OUTPUT_DIR` | Summarizes an existing run |
| `doctor MODEL --write-config-template CONFIG` | Diagnostic/template generation |
| `template --formulation NAME` | Declared contract; alternatively `--material NAME` |
| `replay RECORD --object OBJ --contract MAPPING --out DIR` | Pinned connected J2 replay; `--solve` creates a synthetic equilibrated record; `--verify` adds independent ORIGINAL FD |
| `request --model INP --odb ODB --material OBJ --request JSON --out DIR` | Actual ODB collaborator interface below |

`--odb` on legacy assemble/requirements/verify is an alias for exported
`--fields` JSON, **not** direct binary ODB extraction. Use `request` for the
actual ODB interface. `--no-fd` does not establish independent verification.
Do not use `init --force` on evidence you need to retain.

## Five Reproduced Examples

For a complete repeat with automatic numerical assertions and raw capture:

```sh
"$PY" scripts/audit_recovery_usage.py --umat "$UMAT" --phase examples
```

This runs both repositories' five workflows sequentially in new scratch space.
[Raw results](evidence/usage_examples.json) embed numerical proof and actual
commands; temporary compiler products are not portable inputs. Individual RA
commands below use a new `OUT=$(mktemp -d /tmp/resasm-examples-XXXXXX)`.

### 1. Direct Python Residual

```sh
OUT=$(mktemp -d /tmp/resasm-examples-XXXXXX)
resasm init --template python --out "$OUT/direct"
resasm check "$OUT/direct/resasm.yml"
resasm run "$OUT/direct/resasm.yml"
```

Inputs: [template](../templates/user_python_residual/resasm.yml), cubic spring
R=k*u^3-f, k=2, f=16, u=2. Genuine OTILib order two. Recovered derivatives
du/dk=-1/3, du/df=1/24, d2u/dk2=2/9, d2u/dkdf=-1/144,
d2u/df2=-1/576 all agree with the closed form within `1e-8` absolute.
Raw coefficients are not derivatives: the k-squared coefficient is 1/9.
The captured command used the public `init` on this same template. Missing
OTILib is an error, not a substitute algebra.

### 2. Order-Two Black-Box Residual

```sh
resasm init --template blackbox-order2 --out "$OUT/blackbox"
resasm check "$OUT/blackbox/resasm.yml"
resasm run "$OUT/blackbox/resasm.yml"
```

Inputs: [template](../templates/user_blackbox_order2_residual/resasm.yml), real
executable evaluating R=k^2*u^3-f, k=2, f=32, u=2. The framework receives Taylor
coefficients, solves and applies factorial recovery. Checked analytically:
du/dk=-2/3, du/df=1/48, d2u/dk2=5/9, d2u/dkdf=-1/144,
d2u/df2=-1/2304, absolute tolerance `1e-8`. This executable is a mathematical
example, not a mocked FE solver or proof of another private solver's derivatives.

Both template runs write public summaries/norm CSVs and private direction-map
JSON, RHS and sensitivity NPZ arrays below their `resasm_output` directory.
Private arrays include coefficients, recovered derivatives and recovery factors.
The analytic audit is separate from any pending black-box local FD checks.
GUI equivalent: Run with the copied configuration, through the same CLI bridge;
these two workflows were not independently clicked through in this audit.

### 3. Stress-Driven C3D8

```sh
MODEL=residual_core/examples/minimal_c3d8_stress_driven
resasm assemble "$MODEL/model.json" --mode stress-driven \
  --fields "$MODEL/fields.json" --out "$OUT/stress.npy"
```

One unit cube, eight IPs, supplied S11=100. Measured 24 DOFs,
norm(R)=70.71067811865474 and max(abs(R)) approximately 25, agreeing with analytic face
tractions to `1e-10`. No material tangent or parameter derivative is inferred
from supplied stress. This unloaded example's nonzero residual is expected,
not equilibrium. GUI: Assemble, this model/fields, stress-driven mode.

### 4. Nonlinear J2 C3D8 Replay

```sh
"$PY" scripts/reproduce_imqcam_pipeline.py --skip-abaqus --out "$OUT/j2"
```

Inputs: [cyclic model](../examples/imqcam_j2_cantilever/model.json) and companion
m3_j2 contract. A fresh compiled provider, 7 increments, 8 IPs, four parameters
E,nu,SIGY0,H. The retained manifest says `passed=true`; whole-history
displacement/stress/state derivatives and FD plateau checks pass `2e-6` scaled
error against independently compiled ORIGINAL material and re-equilibrated FE
paths. The record is `synthetic_converged_fe`, not an Abaqus J2 cantilever.
Public plot/CSV and private manifest, provider, record and field arrays remain
under the output. Missing inputs or failed checks exit nonzero. GUI: Advanced
Replay with the generated model/object/contract, Solve and Verify. Existing GUI
integration tests exercise that route; the audit's GUI launch is not a fresh
interactive execution of every example.

### 5. Bounded Finite-Strain C3D8

```sh
"$PY" examples/finite_strain_c3d8/benchmark.py --out "$OUT/finite"
resasm --config "$OUT/finite/config.json" assemble "$OUT/finite/model.json" \
  --mode material-replay --tangent
resasm --config "$OUT/finite/config.json" sensitivity "$OUT/finite/model.json" \
  --params "$OUT/finite/params.json"
```

The benchmark creates two distorted shared-node elements, 36 DOFs, with
mu=2.3, lambda=4.1, manufactured fixed loads and a large rotation. Its current
report passed independent first-Piola force (`1e-12`), complete nodal tangent
FD (`2e-8` at three steps), and genuine OTILib du/dp versus nonlinear re-solves
(`2e-5` at three steps). GUI: Assemble/Sensitivity using the generated config.
The standalone benchmark and both public assembly/sensitivity commands were
executed here and exited zero, in addition to their regression coverage.
Only stateless isotropic total neo-Hookean response is supported. Rotating
plastic history, follower loads, geometric OTI and higher-order finite FE are
not supported. See [full example](../examples/finite_strain_c3d8/README.md).

## Connected Presentation Interface

Developer: build the provider in the companion repository using its existing
`umat-oti-provider build CONTRACT --out NEW_DIR`. Share the compiled object and
completed generated sidecar, not the original source or input transformation
contract. The collaborator supplies exactly these four explicit inputs:

```sh
resasm request --model Analysis.inp --odb Analysis.odb \
  --material OTI_UMAT.obj --request sensitivity_request.json --out new_results
```

Keep the unchanged generated `Mapping.json` beside the object. Discovery uses
the object-stem JSON or Mapping.json; conflicting sidecars require explicit
`--mapping PATH`. Object SHA-256, source fingerprint, dimensions, ABI,
PROPS/direction ordering and derivative/Voigt layout are checked automatically.
The object bundles ORIGINAL and OTI binary routines. Source is not required;
the generated link shim is not the private constitutive source. Trust the binary
supplier: source-denial audit hooks are not an OS sandbox or binary certification.

Exactly three public files are written at the output root:

| File | Meaning |
| --- | --- |
| `sensitivity_results.json` | Requested scalar values and derivatives, resolved scope and execution metadata |
| `sensitivity_tables.csv` | One row per output/increment/parameter; four outputs x four parameters gives 16 rows here |
| `run_report.txt` | Execution checks, tolerances, verification status and limitations |

Full fields, K/R/derivatives, generated replay record, exporter log and linked
library stay in `private/`. Do not publish that directory as a scalar report.
No transform, source read or production solve occurs during consumption.
`--abaqus PATH` selects licensed odbAccess extraction; missing executable or
failed export is an explicit error.

Request keys are `outputs`, `parameters`, `domain`, `increments`. For example:

```json
{"outputs":[{"name":"loaded_U1","field":"U","component":1,"reduction":"mean"}],"parameters":["E","SIGY0","H"],"domain":{"nodes":[2,3,6,7]},"increments":"LAST"}
```

Fields U/RF use nodes; S/SDV use elements and all eight IPs. Components are
one-based or ALL; stress order is 11,22,33,12,13,23; SDV1 is equivalent plastic
strain. Parameters are a unique subset of E,nu,SIGY0,H. Increments are ALL,
LAST or a unique one-based list; all preceding history is replayed first.
Reductions: component (one location), unweighted sum/mean, Euclidean L2, signed
max. Zero-norm L2 and tied maxima fail. Mean is not volume weighted, sum is not
a volume integral, and L2 is not von Mises. Unknown fields/ids/options fail.

Supported INP/ODB: one homogeneous pinned J2 C3D8/B-bar static NLGEOM=NO step,
one untransformed instance, virgin state, zero fixed BCs, ramped nodal loads.
Every frame including frame zero must contain U/RF/CF/S/SDV1 and matching
mesh/history. No interpolation of missing frames. Pressure, body loads, contact,
amplitudes, multiple steps/materials/instances, nonzero prescribed displacement,
initial state, finite strain, generic FCC and full-size cantilevers are refused.

ODB acceptance: scaled free residual <`1e-5`; stress/RF relative tolerance
`2e-5` with field-scaled absolute floor; state `rtol=2e-5, atol=1e-8`; zero BC
residues below `1e-12*mesh_extent` only. Float32 ODB precision limits accuracy.
Ordinary execution reports `verified=false`; `--validate` explicitly invokes
independent ORIGINAL FD and may fail its stricter double-precision gate.

Fresh source-denied consumption passed the separate uniaxial J2 analytic check
(`2e-5` relative for nonzero derivatives; `1e-8` absolute for zero nu derivative).
[Raw connected evidence](evidence/usage_presentation.json) embeds input hashes,
the exact public file list and the independent check. Reproduce with:

```sh
"$PY" scripts/audit_recovery_usage.py --umat "$UMAT" --phase presentation
```

This copies only the five existing genuine collaborator artifacts into new
scratch space and calls `consume`, never `prepare`. The external ODB is required
and is not shipped in the repository. See [exact interface](PRESENTATION_INTERFACE.md).

## GUI And Report Interpretation

```sh
"$PY" -m streamlit run scripts/app.py --server.address=127.0.0.1 --server.port=8501
```

Choose a free port and open the printed URL. Sensitivity Request is the primary
screen with the four input paths/uploads, optional advanced mapping and three
downloads. Legacy workflows remain available through the other screens. For an
installed wheel, locate `residual_core.app.streamlit_app` with importlib.util
and pass its file path to Streamlit, as in the installation evidence.

[GUI evidence](evidence/usage_gui.json): primary render without exception and
actual HTTP-ready server, stopped by owned PID. No audit URL remains live.
Existing real-browser evidence is available as local artifacts only:
[report](../.pytest_cache/presentation_browser_cli_final/browser_report.json),
[desktop](../.pytest_cache/presentation_browser_cli_final/desktop.png),
[mobile](../.pytest_cache/presentation_browser_cli_final/mobile.png).
Those local cache artifacts are not promised in a clean clone; no large image
was copied into this report. The recorded browser run checked actual ODB
execution and all three downloads, not just initial rendering.

Read verdicts separately: command executed; residual assembled; free equilibrium
checked/passed; tangent available/independently verified; derivative calculated/
independently verified; reference resolved; Abaqus comparison available. A full
residual includes reactions and need not be zero. An available tangent, finite
output or successful sensitivity linear solve is not independent verification.
Historical fixture checks can say held, failed or not established; exact zeros
and unresolved comparisons are not automatically independent agreements.

## Requirement Coverage

### Review Follow-Up: Failure Privacy

Review checks on 2026-09-18: **60 passed, zero failures/skips**, comprising 52
presentation tests and eight unchanged thin-CLI GUI guards. Retained result:
[focused JUnit](evidence/review_fixes_focused.xml). From the workspace root,
with the recovery import and OTILib environment above:

```sh
.venv/bin/pytest -q -ra imq-ra-recovery/tests/integration/test_presentation_request.py imq-ra-recovery/tests/framework/test_gui_is_a_thin_cli_front_end.py --junitxml=imq-ra-recovery/docs/evidence/review_fixes_focused.xml
```

Failure output keeps the exact public name `run_report.txt`. It reports a fixed
category/action and the relative diagnostic path `private/error_report.txt`;
the latter retains the original exception and traceback. CLI failures return 2
and the thin GUI displays the sanitized CLI diagnostic. Missing input paths are
identified by role without echoing private path text; missing ODB field actions
name only allowlisted fields. If the output directory cannot be prepared or
written, the CLI reports that private diagnostics are unavailable.

Two explicitly marked unit cases inject solver failure and exporter subprocess
stdout/stderr with a recognizable sensitive marker, checking service errors,
all public files that exist, CLI stderr/stdout and GUI-visible text. They are
disclosure tests, **not numerical integration proof**. Four missing-path unit
cases and five missing-field transport cases also pass. The compiled numerical
presentation regressions remain unchanged, including the genuine archived J2
field check; no core mathematics was mocked in those checks. No full suite,
fresh licensed analysis or clean-install gate was run for this review.

The [274-row ledger](COMPLETION_LEDGER.md) and [machine-readable requirement
index](evidence/usage_requirements.json) distinguish 104 implemented bounded
rows, 159 partial-evidence rows and 11 unestablished release-gate rows.
**Zero rows have final clean-install PASS; all 274 remain open
under that rule.** No source-suite or working-tree wheel result is promoted.
Highest priorities: final committed-pair clean-clone reproduction; all-example
GUI execution and installed example discovery; general stateful finite/FCC and
full-size prescribed-displacement histories; higher-order full FE; current
corpus evidence and missing historical cross-reader checks. A diagnostic refusal
is safer than silent physics loss, but it does not implement the refused feature.