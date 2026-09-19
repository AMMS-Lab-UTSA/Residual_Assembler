# Connected workflow: from UMAT source to sensitivities

This page describes the developer-side connected workflow: UMAT-OTI transforms
and compiles a UMAT into a provider, and `resasm replay` links that provider to
solve for total-history sensitivities of a bounded J2 model, with a
one-command reproducer that checks the whole chain against independent
references. It is for material developers and reviewers who want to see the
provider-to-sensitivity pipeline end to end.

Collaborators who hold only the four shared files (`Analysis.inp`,
`Analysis.odb`, `OTI_UMAT.obj` with `Mapping.json`, `sensitivity_request.json`)
do not need this page: they run `resasm request`
([REQUEST_INTERFACE.md](REQUEST_INTERFACE.md)), which exports every ODB
increment with licensed Abaqus Python, checks that the deck and the ODB agree,
and hands full-size models to `resasm history`
([REPLAY_HISTORY.md](REPLAY_HISTORY.md)).

## What the reproducer runs

This is a runnable integration of the two products on a bounded model, not
the full-size cantilever (for that, see
[examples/cantilevers](../examples/cantilevers/README.md)). UMAT-OTI freshly
transforms and compiles the ORIGINAL `m3_j2` UMAT. Residual_Assembler links
that object, loads a C3D8 model and record, marches eight integration points
through seven increments, assembles `K`, `R` and `Rp`, and solves the
constrained total-history sensitivity problem for E, nu, SIGY0 and H.

The J2 example is a **synthetic converged small-strain FE calculation**, not an
Abaqus export. A separate, genuine Abaqus 2021 export of a non-uniform
elastic C3D8 model (`tests/abaqus_derivative_export/`) independently checks
B-bar assembly, integration-point order, equilibrium and sensitivities. No
prebuilt provider or generated scratch directory is a required input.

## Environment

Run from a Residual_Assembler checkout with the UMAT_source_transformation
checkout beside it and a genuine OTILib build (`scripts/setup_otilib.sh`, or
[OTILIB_VENV.md](OTILIB_VENV.md)):

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

The complete reproduction and its tests need Python 3.11, gfortran, NumPy,
matplotlib, pytest and Streamlit. The compiled provider uses the
transformer's bundled Fortran OTI implementation; the separate OTILib build is
exercised by the Python scalar guards and the OTILib tests. Neither Dual1 nor
the unrelated PyPI `pyoti` is a substitute.

## One command

Choose an output directory that does not exist:

```sh
"$PY" scripts/reproduce_connected_pipeline.py \
  --provider-repo "$UMAT_OTI_REPO" \
  --skip-abaqus --out .pytest_cache/connected_reproduction
```

The script transforms and compiles anew, invokes the public CLI twice (solve,
then replay), checks the ORIGINAL material paths and independently
equilibrated ORIGINAL models at three centred finite-difference steps, and
verifies the genuine elastic Abaqus fixture. It produces:

- `manifest.json`: success or failure and output locations; a failure exits 2.
- `private/manifest.json`: both Git commits, dirty status, per-source SHA-256
  and source-tree hashes, interpreter, environment and dependencies, compiler,
  original, object and contract hashes, exact subprocess commands and the
  numerical checks.
- `private/provider/`: the new object, contract and generated Fortran and
  build sources.
- `private/solve/private/`: the converged record and the full `K`, `R`, `Rp`,
  state and derivative arrays.
- `private/replay/`: the independently loaded replay results.
- `public/verification.json`, `public/abaqus_fixture.json`,
  `public/fd_table.csv` and `public/j2_history.png`: portable summaries, a table
  and a plot.

`--skip-abaqus` is allowed because the elastic exported fields are committed
in the repository, not synthesized for this run. Without that flag the script
refuses explicitly: it does not launch licensed J2 jobs. The flag never
changes the synthetic provenance of the J2 result. The availability of
`abaqus` is recorded; availability alone is not evidence of a licensed run.

## Public CLI

```sh
"$PY" -m umat_oti.provider build \
  "$UMAT_OTI_REPO/parameter_sensitivity/models/m3_j2/contract_v2.json" \
  --out .pytest_cache/my_j2_provider

"$PY" -m residual_core.ui.cli replay \
  examples/bounded_j2_c3d8/model.json \
  --object .pytest_cache/my_j2_provider/umat_m3_j2_oti.obj \
  --contract .pytest_cache/my_j2_provider/umat_m3_j2_oti.json \
  --solve --verify --out .pytest_cache/my_j2_solve

"$PY" -m residual_core.ui.cli replay \
  .pytest_cache/my_j2_solve/private/record.json \
  --object .pytest_cache/my_j2_provider/umat_m3_j2_oti.obj \
  --contract .pytest_cache/my_j2_provider/umat_m3_j2_oti.json \
  --verify --out .pytest_cache/my_j2_replay
```

`--solve` generates a synthetic equilibrated record. Without it, every
increment must contain a converged `u`; the CLI checks equilibrium, prescribed
DOFs, and recorded stress and state when supplied. Omit `--verify` only to
execute without finite differences; the summary then says `verified=false`.
Invalid input, failed convergence, a fingerprint mismatch or failed
verification returns 2 and writes a failed public summary. Do not treat old
private output from a failed rerun as valid.

The **Advanced Replay** GUI tab is a thin call to this same command, using the
GUI's usual exit-code and reporting machinery. Launch it in the same
environment:

```sh
"$PY" -m streamlit run scripts/app.py
```

Enter the new object and contract paths and choose Solve model and Verify. The
GUI regression test clicks this button and runs the compiled material; it does
not mock the CLI or the constitutive mathematics.

## Sensitivity semantics

`PathMaterial.march_oti` and `march_fast` compute parameter derivatives along
a **fixed strain history**, carrying stress and physical state derivatives. A
final solve with only those derivatives is not, in general, the derivative of
a path-dependent equilibrium solution.

For each converged increment `n`, the connected replay forms the current-`u`
partial by carrying the preceding total stress and state derivatives and
seeding `d(delta_strain)/dp = -B du[n-1]/dp`. It assembles `Rp_history` and
solves `K_ff du[n]/dp = -Rp_history_f`. It then recomputes the update
derivative with `B (du[n]/dp - du[n-1]/dp)` and commits the resulting total
stress and state sensitivities. Physical history is never committed during
Newton trials.

The pinned J2 update depends on strain through
`stress_trial = stress_previous + C_elastic delta_strain`. Seeding
`dstress_previous + C_elastic d(delta_strain)/dp` in the 18-argument EVAL ABI
is therefore exactly the strain-history chain rule, including the state
derivative. The ORIGINAL UMAT supplies `C_elastic` at virgin zero strain;
Residual_Assembler contains no replacement J2 return mapping. The provider
adds the direct parameter terms, including `dC_elastic/dp * delta_strain`.

This equivalence is **not a property of the generic UMAT ABI**. The command
checks the known `m3_j2` source fingerprint, dimensions, parameter order,
EVAL/MARCH argument counts, linked MARCH availability and object hash, and
requires matching source provenance in the record. Other materials are
rejected here; the history engine handles any provider through the
`UMAT_OTI_EVAL_TOTAL` entry point instead ([REPLAY_HISTORY.md](REPLAY_HISTORY.md)).

The outputs `du_dp`, `dsigma_dp` and `dstate_dp` are total
equilibrated-history values. `fixed_path_dsigma_dp`, `fixed_path_dstate_dp`,
`Rp_fixed_path` and `fixed_path_local_solve` are separately labelled
diagnostics; the last is not presented as a solution derivative.
`Rp_history` holds the current displacement fixed but includes the dependence
on previous solutions. `K` is the current increment's algorithmic tangent. `R`
includes reactions on prescribed DOFs.

## Verification and limits

The finite-difference reference calls the bundled untransformed ORIGINAL UMAT,
not the lifted update. Each perturbed FE model is restarted from the virgin
state and independently converged with ORIGINAL stresses and tangents.
Fixed-path finite differences independently march the ORIGINAL at every
integration point. All four parameters, seven increments and eight
integration points are checked for displacement, stress and state, with
branch-preserving perturbations `h = (1e-4, 3e-5, 1e-5) * max(abs(p),1)`. Every
step must meet a scaled error of `2e-6`, and the last two estimates must also
meet that plateau limit. The scale is the largest reference magnitude over the
complete array for that parameter, floored at `1e-12`; it is not a relative
error of individual near-zero entries.

Newton uses a residual-decreasing line search and a normalized free-residual
limit of `1e-14`. A singular tangent, a failed line search or the 40-iteration
limit fails explicitly. This is a dense, bounded demonstration, not a scalable
solver.

`resasm replay` supports: homogeneous pinned J2, first derivatives,
small-strain C3D8 with full or B-bar integration, virgin initial state, zero
prescribed displacements, and parameter-independent nodal loads, geometry and
boundaries. Unsupported context and record fields fail rather than
disappearing from equilibrium. The legacy shim does not carry temperature,
coordinates, rotation, time-dependent constitutive context or KINC; the pinned
source does not use them. Live OTI values are rejected at real-field and ABI
boundaries rather than silently cast. ABI doubles that have already been
extracted are valid inputs to NumPy assembly and linear solves.

This command makes no claim about finite strain, shape, load or boundary
derivatives, multiphysics, arbitrary UMATs, nonzero initial histories, higher
derivatives or cross-platform binaries; full-size models and other providers
are the history engine's job. The deliberately strong cyclic example can
accumulate about 10% displacement and strain on a unit cube: it validates the
specified small-strain algebra, not its physical accuracy at large
deformation. This synthetic reproduction makes no J2 Abaqus comparison;
`resasm request` provides the genuine J2 ODB comparison
([REQUEST_INTERFACE.md](REQUEST_INTERFACE.md)).

The exact commands, measured errors, hashes and remaining dataset skips are in
the [integration evidence](evidence/recovery_integration.md).
