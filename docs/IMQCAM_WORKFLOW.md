# Presentation-Matching Collaborator Workflow

The primary interface now consumes the original slide-12 inputs directly:

```sh
cd /path/to/Residual_Assembler
export PY=/path/to/venv/bin/python        # a venv with both packages installed
"$PY" -m residual_core.ui.cli request \
  --model /path/to/Analysis.inp --odb /path/to/Analysis.odb \
  --material /path/to/OTI_UMAT.obj \
  --request /path/to/sensitivity_request.json --out /path/to/new_results
```

It discovers the generated mapping, exports every ODB increment with licensed
Abaqus Python, checks INP/ODB agreement, and uses the total-history J2 engine.
No source transformation, synthetic solution, production rerun or manually
created record is required. Outputs are exactly `sensitivity_results.json`,
`sensitivity_tables.csv`, `run_report.txt`; full fields/records remain private.
FD is off by default. The primary GUI **Sensitivity Request** screen has the
four selectors, Advanced mapping selection, and exact output downloads.

Scope: one small-strain static C3D8 step, homogeneous pinned m3_j2, virgin state,
zero BCs, ramped nodal loads, complete U/RF/CF/S/SDV1 history. A genuine
four-increment Abaqus J2 job and source-denied collaborator run are now measured.
See [PRESENTATION_INTERFACE.md](PRESENTATION_INTERFACE.md) for format and limits;
[evidence/recovery_presentation.md](evidence/recovery_presentation.md) for proof.

## Advanced Developer And Replay Workflow

The earlier workflow below is retained for developer validation, not ordinary
collaborator use. Its synthetic J2 fixture and archived elastic comparison are
separate from the new genuine J2 ODB proof above.

This is a runnable Program 1 -> Program 2 integration, not the full IMQCAM
cantilever claim. Program 1 freshly transforms and compiles the ORIGINAL
`m3_j2` UMAT. Program 2 links that object, loads a C3D8 model/record, marches
eight integration points through seven increments, assembles K/R/Rp, and solves
the constrained total-history sensitivity problem for E, nu, SIGY0 and H.

The shipped J2 example is a **synthetic converged small-strain FE calculation**,
not an Abaqus export. A separate genuine, archived Abaqus 2021 elastic fixture
independently validates B-bar assembly, IP order, equilibrium and sensitivities.
No prebuilt provider or generated scratch directory is a required input.

## Environment

Run from a Residual_Assembler checkout with the UMAT_source_transformation
checkout beside it and a genuine OTILib build (`scripts/setup_otilib.sh`):

```sh
cd /path/to/Residual_Assembler
export PY=/path/to/venv/bin/python
export UMAT_OTI_REPO=/path/to/UMAT_source_transformation
export OTILIB=/path/to/otilib/build_py311
export PYTHONPATH="$PWD:$UMAT_OTI_REPO/src:$OTILIB"
export PYOTI_PATH="$OTILIB"
export OTILIB_ROOT="$OTILIB"
export RUN_OTILIB_TESTS=1
"$PY" -c 'import residual_core,umat_oti,pyoti; print(residual_core.__file__); print(umat_oti.__file__); print(pyoti.__file__)'
```

Requires the W3 recovery provider sources, Python 3.11, gfortran, NumPy,
matplotlib, pytest and Streamlit for the complete reproduction/test workflow.
The compiled provider uses the transformer's bundled Fortran OTI implementation;
the separate real OTILib build is exercised by the Python scalar guards and
existing W1 tests. It is not substituted by Dual1 or the unrelated PyPI `pyoti`.

## One Command

Choose an output directory that does not exist:

```sh
"$PY" scripts/reproduce_imqcam_pipeline.py \
  --provider-repo "$UMAT_OTI_REPO" \
  --skip-abaqus --out .pytest_cache/imqcam_reproduction
```

The script transforms and compiles anew, invokes the public CLI twice (solve
then replay), checks ORIGINAL material paths and independently equilibrated
ORIGINAL models at three centered FD steps, and verifies the real elastic
Abaqus fixture. It produces:

- `manifest.json`: success/failure and output locations; failure exits 2.
- `private/manifest.json`: both Git commits, dirty status, per-source SHA-256
  hashes and source-tree hashes, interpreter/environment/dependencies, compiler,
  original/object/contract hashes, exact subprocess commands and numerical checks.
- `private/provider/`: new object, contract and generated Fortran/build sources.
- `private/solve/private/`: converged record and full K/R/Rp, state and derivative arrays.
- `private/replay/`: independently loaded replay results.
- `public/verification.json`, `public/abaqus_fixture.json`, `public/fd_table.csv`
  and `public/j2_history.png`: portable summaries, table and plot.

`--skip-abaqus` is allowed because the elastic exported fields were already
committed in upstream history, not synthesized for this run. Without that flag
the bounded script refuses explicitly: fresh licensed J2 jobs are not implemented.
The flag never changes the J2 result's synthetic provenance. Availability of
`abaqus` is recorded; availability alone is not evidence of a licensed run.

## Public CLI

```sh
"$PY" -m umat_oti.provider build \
  "$UMAT_OTI_REPO/parameter_sensitivity/models/m3_j2/contract_v2.json" \
  --out .pytest_cache/my_j2_provider

"$PY" -m residual_core.ui.cli replay \
  examples/imqcam_j2_cantilever/model.json \
  --object .pytest_cache/my_j2_provider/umat_m3_j2_oti.obj \
  --contract .pytest_cache/my_j2_provider/umat_m3_j2_oti.json \
  --solve --verify --out .pytest_cache/my_j2_solve

"$PY" -m residual_core.ui.cli replay \
  .pytest_cache/my_j2_solve/private/record.json \
  --object .pytest_cache/my_j2_provider/umat_m3_j2_oti.obj \
  --contract .pytest_cache/my_j2_provider/umat_m3_j2_oti.json \
  --verify --out .pytest_cache/my_j2_replay
```

`--solve` generates a synthetic equilibrated record. Without it every increment
must contain converged `u`; the CLI checks equilibrium, prescribed DOFs, and
recorded stress/state when supplied. Omit `--verify` only for execution without
FD; the summary then explicitly says `verified=false`. Invalid input, failed
convergence, fingerprint mismatch or failed verification returns 2 and writes a
failed public summary. Do not treat old private output from a failed rerun as valid.

The **Advanced Replay** GUI tab is a thin call to this same command, using the existing
GUI exit-code/reporting machinery. Launch manually with the same environment:

```sh
"$PY" -m streamlit run scripts/app.py
```

Enter the new object/contract paths and choose Solve model and Verify. The GUI
regression actually clicks this button and runs the compiled material; it does
not mock the CLI or constitutive mathematics. No background GUI server was
started during the synchronous recovery invocation.

## Sensitivity Semantics

`PathMaterial.march_oti`/`march_fast` compute parameter derivatives along a
**fixed strain history**, carrying stress and physical state derivatives. A
final solve with just those derivatives is not generally the derivative of a
path-dependent equilibrium solution.

For each converged increment n, connected replay forms the current-u partial
by carrying the preceding total stress/state derivatives and seeding
`d(delta_strain)/dp = -B du[n-1]/dp`. It assembles `Rp_history` and solves
`K_ff du[n]/dp = -Rp_history_f`. It then recomputes the update derivative using
`B (du[n]/dp - du[n-1]/dp)` and commits the resulting total stress/state
sensitivities. Physical history is never committed during Newton trials.

The pinned J2 update depends on strain through
`stress_trial = stress_previous + C_elastic delta_strain`. Therefore seeding
`dstress_previous + C_elastic d(delta_strain)/dp` in the 18-argument EVAL ABI
is exactly the strain-history chain rule, including the state derivative.
The ORIGINAL UMAT supplies `C_elastic` at virgin zero strain; Program 2 contains
no replacement J2 return mapping. The provider automatically adds the direct
parameter terms, including `dC_elastic/dp * delta_strain`.

This equivalence is **not a generic UMAT ABI property**. The command checks the
known m3_j2 source fingerprint, dimensions, parameter order, EVAL/MARCH argument
counts, linked MARCH availability and object hash, and requires matching source
provenance in the record. Other materials are rejected.

Output `du_dp`, `dsigma_dp`, `dstate_dp` are total equilibrated-history values.
`fixed_path_dsigma_dp`, `fixed_path_dstate_dp`, `Rp_fixed_path` and
`fixed_path_local_solve` are separately labelled diagnostics; the latter is not
advertised as a solution derivative. `Rp_history` holds current displacement
fixed but includes previous solution dependence. K is the current increment's
algorithmic tangent. R includes reactions on prescribed DOFs.

## Verification and Limits

The FD reference calls the bundled untransformed ORIGINAL UMAT, not the lifted
update. Each perturbed FE model is restarted from virgin state and independently
converged with ORIGINAL stresses/tangents. Fixed-path FD independently marches
the ORIGINAL at every IP. All four parameters, seven increments and eight IPs
are checked for displacement/stress/state, with branch-preserving perturbations
`h = (1e-4, 3e-5, 1e-5) * max(abs(p),1)`. Every step must meet scaled error
`2e-6`; the last two estimates must also meet that plateau limit. Scaling is the
maximum reference magnitude over the complete array for that parameter, floored
at `1e-12`, not a relative error for individual near-zero entries.

Newton uses a residual-decreasing line search and a normalized free residual
limit of `1e-14`. A singular tangent, failed line search or 40-iteration limit
fails explicitly. This is a dense, bounded demonstration, not a scalable solver.

Supported: homogeneous pinned J2, first derivatives, small-strain C3D8 with
full or B-bar integration, virgin initial state, zero prescribed displacements,
parameter-independent nodal loads/geometry/boundaries. Unsupported context and
record fields fail rather than disappearing from equilibrium. The legacy shim
does not carry temperature, coordinates, rotation, time-dependent constitutive
context or KINC; the pinned source does not use them. Live OTI values are rejected
at real-field/ABI boundaries rather than silently cast. Already extracted ABI
doubles are valid inputs to NumPy assembly and linear solves.

No finite-strain, shape/load/BC derivatives, multiphysics, arbitrary UMATs,
nonzero initial histories, higher derivatives, FCC/full-size cantilevers,
cross-platform binaries or clean-install claims are made. The deliberately
strong cyclic example can accumulate about 10% displacement/strain on a unit
cube; it validates the specified small-strain algebra, not its physical accuracy
at large deformation. No J2 Abaqus comparison was performed by this older
synthetic reproduction. The presentation workflow above now provides a genuine
bounded J2 ODB comparison, not a full-size cantilever claim.

See [evidence/recovery_integration.md](evidence/recovery_integration.md) for exact
commands, measured errors, hashes and remaining full-suite dataset skips.