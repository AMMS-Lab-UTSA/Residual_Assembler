# Example 3: the four-file request on one element

The complete collaborator workflow on the smallest possible Abaqus analysis:
one C3D8 element of a J2-plastic material, pulled into the plastic range.
From `Analysis.inp`, `Analysis.odb`, the compiled material (`OTI_UMAT.obj`
with its `Mapping.json`) and `sensitivity_request.json`, one command returns
the derivatives of displacement, reaction, stress and plastic strain with
respect to all four material constants. Every number has a closed-form
reference.

| | |
| --- | --- |
| Folder | `examples/presentation_request/` (the deck and the request) |
| Needs Abaqus | Yes: once to create the ODB, and Abaqus Python to read it |
| Needs OTILib | No |
| Run time | about 1 s (plus about 6 s to build the material provider once) |

## What it demonstrates

This is the interface meant for the person who ran the analysis. They do not
need the material's source code and they do not re-run the analysis. They
hand four files to `resasm request` and receive three public files:
`sensitivity_results.json`, `sensitivity_tables.csv` and `run_report.txt`.

The model is small enough that the answer is known exactly, so the example is
also a check of the whole chain: the ODB export, the replay of the compiled
material at all eight integration points and four increments, the assembly of
the residual and its tangent, and the sensitivity solve.

## The mathematics

**The replay.** At every increment `n` the residual is assembled from the
recorded displacements and the replayed material:

    R_n(u_n, p) = sum_e sum_q w_q B_q^T sigma_{n,q}(u, history, p) - F(t_n).

Differentiating `R_n = 0` with respect to the material parameters `p` gives
the sensitivity equation solved at every increment, on the free degrees of
freedom:

    K du_n/dp = - dR_n/dp,     K = dR_n/du  (assembled from the material tangent DDSDDE).

The right-hand side carries the whole history: the stress and state
derivatives of increment `n-1` are passed into the material at increment `n`.
The compiled material returns the derivatives exactly, through OTI numbers.

**The reference.** The element is the unit cube, supported on the planes
`x = 0`, `y = 0` and `z = 0`, and loaded by 75 per node on the four nodes of
the face `x = 1` (node set `LOADED`), ramped over four increments. The state is
uniaxial stress. At the last increment `sigma = 4 x 75 / 1 = 300`, above the
initial yield stress `SIGY0 = 250`, and with linear hardening `H`:

    eqplas = (sigma - SIGY0) / H = 0.025,        U1 = sigma / E + eqplas = 0.0264285714...
    dU1/dE = -sigma / E^2 = -6.8027e-09,         dU1/dnu = 0,
    dU1/dSIGY0 = -1 / H = -5.0e-04,              dU1/dH = -(sigma - SIGY0) / H^2 = -1.25e-05.

The loading is force-controlled and the structure is statically determinate,
so the stress (300) and the support reaction (-300) do not depend on the
material at all: their derivatives are zero. The plastic strain has the same
`SIGY0` and `H` derivatives as `U1` and none with respect to `E` and `nu`.
Material constants: `E = 210000`, `nu = 0.3`, `SIGY0 = 250`, `H = 2000`.

## Inputs

| File | Where it comes from |
| --- | --- |
| `Analysis.inp` | this folder. One C3D8, one `*User Material` with 4 constants and 1 state variable, three symmetry supports, a `*Cload`, one static step of four increments, field output `U, RF, CF, S, SDV` every increment |
| `sensitivity_request.json` | this folder. Four outputs at the last increment: the mean `U1` of node set `LOADED` (nodes 2, 3, 6, 7), the summed `RF1` of the support nodes 1, 4, 5, 8, the mean `S11` and the mean `SDV1` (equivalent plastic strain) over the element; all four parameters |
| `Analysis.odb` | not shipped (a binary file). Produced once by Abaqus from `Analysis.inp` with the ORIGINAL J2 UMAT, see step 2 |
| `OTI_UMAT.obj`, `Mapping.json` | built from the companion repository with `umat-oti-provider build`, see step 1 |

## Run it from the command line

From the Residual_Assembler root, with both repositories side by side and the
environment of [INSTALL.md](../../docs/INSTALL.md):

```bash
export RA="$PWD"
export UMAT="$RA/../UMAT_source_transformation"
export WORK="$HOME/resasm_work"; mkdir -p "$WORK"
```

**Step 1. Build the material provider** (the material owner does this once):

```bash
umat-oti-provider build "$UMAT/parameter_sensitivity/models/m3_j2/contract_v2.json" \
    --out "$WORK/provider_j2"
```

It prints the paths of `umat_m3_j2_oti.obj` (the compiled material, both the
ORIGINAL routine and its OTI version) and `umat_m3_j2_oti.json` (the completed
mapping: parameter names, PROPS positions, layouts and the SHA-256 of the
object). It also leaves a `build-*` folder with the generated sources; that
folder stays with the material owner. Only the `.obj` and the `.json` are
handed over.

**Step 2. Put the four inputs in one folder.**

```bash
D="$WORK/one_element"; mkdir -p "$D"
cp examples/presentation_request/Analysis.inp examples/presentation_request/sensitivity_request.json "$D/"
cp "$WORK/provider_j2/umat_m3_j2_oti.obj"  "$D/OTI_UMAT.obj"
cp "$WORK/provider_j2/umat_m3_j2_oti.json" "$D/Mapping.json"
```

If you already have the ODB of this deck, copy it to `"$D/Analysis.odb"`.
Otherwise create it once with Abaqus and the ORIGINAL UMAT (a Fortran compiler
configured for Abaqus user subroutines is needed):

```bash
cd "$D"
abaqus job=Analysis input=Analysis.inp \
    user="$UMAT/parameter_sensitivity/models/m3_j2/umat.for" cpus=1 interactive
grep "COMPLETED SUCCESSFULLY" Analysis.sta
```

This Abaqus command was not re-run for this guide. The measured results below
use an ODB made with Abaqus 2021.HF5 from this deck and this UMAT
(SHA-256 `54745b97b515de159c3251e8a3ba0ed89ea1c65452e388ef34c8c81794e3b736`).

**Step 3. Run the request.**

```bash
cd "$D"
resasm request --model Analysis.inp --odb Analysis.odb \
    --material OTI_UMAT.obj --request sensitivity_request.json --out results
```

Measured output (exit 0):

```text
request executed: 4 scalar results; verified=False
```

`verified=False` is correct and expected: an ordinary run computes the
derivatives but does not check them against an independent method. The
checks are described below.

## Run it in the GUI

Start the GUI with `streamlit run scripts/app.py` from the repository root
([GUI guide](../../docs/GUI_GUIDE.md)). It opens on **Sensitivity Request**.

1. Fill the three path fields on the right: **OTI_UMAT.obj path**,
   **Analysis.inp path** and **Analysis.odb path** (the files in `"$D"`), or use
   the **Upload** buttons on the left. Leave **sensitivity_request.json**
   empty to build the request with the controls below, or give the file to use
   it as written.
2. Under **3. Outputs and parameters**, the tick boxes `E`, `nu`, `SIGY0`, `H`
   appear, read from the `Mapping.json` beside the object. Leave all four
   ticked.
3. Choose **output** = `displacement U1`, **region** = `node set LOADED`,
   **summary over the region** = `mean`, **increments** = `last increment`.
4. Set **Output directory** to a new folder, for example `$WORK/one_element/gui_results`.
5. Press **Solve**.

Measured on 2026-09-18 (Solve took 1.0 s): the screen reported
`Executed: 1 scalar results. Independent validation: not run.`, offered the
three downloads `sensitivity_results.json`, `sensitivity_tables.csv` and
`run_report.txt`, and showed the full field of the region, one row per node:

| location | value | d/dE | d/dnu | d/dSIGY0 | d/dH | governing parameter |
| --- | --- | --- | --- | --- | --- | --- |
| node 2 | 0.026429 | -6.802720e-09 | 1.172950e-09 | -5.000000e-04 | -1.250000e-05 | SIGY0 |
| nodes 3, 6, 7 | the same values | | | | | SIGY0 |

with **stress reproduced vs. the .odb** = `4.1e-07` and **solution
sensitivity ‖du/dp‖** = `0.00422`. The screenshot shows this run:

![Sensitivity Request after Solve, U1 on node set LOADED](../../docs/screenshots/resasm_solve.png)

Choosing **output** = `state variable SDV1` and **region** = `whole mesh`
gives one row per integration point, eight in all, each with
`EQPLAS = 0.025`, `d/dSIGY0 = -5.000000e-04` and `d/dH = -1.250000e-05`
([screenshot](../../docs/screenshots/resasm_solve_sdv_all_points.png)).
The request actually used is stored in `sensitivity_results.json` under
`request`, so any GUI run can be repeated on the command line.

## What it writes

| Path in `results/` | Contents | Share? |
| --- | --- | --- |
| `sensitivity_results.json` | the request as used, the resolved scope (4 increments replayed, output increment 4, C3D8, selective-reduced integration), metadata (SHA-256 of all five inputs, checks, ODB tolerances, limitations) and the four results with value and derivatives | yes |
| `sensitivity_tables.csv` | one row per output, increment and parameter: 16 rows here | yes |
| `run_report.txt` | what was executed and checked, tolerances, limits | yes |
| `private/fields.json` | the ODB export (every frame: U, RF, CF, S, SDV) | no |
| `private/export_command.json`, `private/odb_export.log` | the Abaqus Python command and its log | no |
| `private/result.json` | full fields and their derivatives at every increment (`u`, `du_dp`, `stress`, `dsigma_dp`, `state`, `dstate_dp`) | no |
| `private/link/` | the small generated link library for the object | no |

The three public files contain no full arrays. `private/` contains the model
and should not be sent with the results.

## Expected output (measured on 2026-09-18)

`run_report.txt`:

```text
Status: executed successfully
Sensitivity semantics: total equilibrated history
Production analysis rerun: no
Material source read/transformed: no
Independent derivative verification: NOT RUN (execution is not independent verification)
History increments replayed: 4
Maximum scaled free residual: 4.12005847e-07 (limit 1e-5)
Checks: mesh ids/coordinates/connectivity; complete virgin-to-final history; time/ramped CF loads; prescribed displacements; free equilibrium; IP stress and state; reactions
Solver: constrained dense tangent solve; reactions Rp_history + K du/dp
```

The results, against the closed form:

| Output | Value | d/dE | d/dnu | d/dSIGY0 | d/dH |
| --- | --- | --- | --- | --- | --- |
| mean_loaded_U1 | 0.026428570970892906 | -6.802719961317867e-09 | 1.1729503437354757e-09 | -4.999999999999996e-04 | -1.249999988950807e-05 |
| reference | 0.02642857142857143 | -6.802721088435374e-09 | 0 | -5.0e-04 | -1.25e-05 |
| support_RF1 | -299.9998763982459 | 1.1e-19 | -4.3e-14 | 3.8e-15 | 1.3e-16 |
| mean_S11 | 299.9998763982459 | -9.2e-19 | 8.3e-15 | 8.3e-15 | 1.4e-16 |
| mean_EQPLAS | 0.024999999779016152 | 2.5e-23 | -7.8e-18 | -4.999999999999996e-04 | -1.2499999889508068e-05 |

## How the result is checked independently

1. **Closed form.** Relative errors of the nonzero `U1` derivatives:
   `E` 1.7e-7, `SIGY0` 8.7e-16, `H` 8.8e-9. The derivatives that should be
   zero are at most 4.3e-14 for the reaction and the stress.
2. **Where the small differences come from.** The ODB stores single-precision
   numbers, and `resasm request` replays the recorded displacements without
   re-equilibrating them (maximum scaled free residual 4.1e-7). That is the
   source of `dU1/dnu = 1.17e-9` (1.3e-8 of `U1` after scaling by `nu`) and of
   the 1.7e-7 in `dU1/dE`. To remove it, replay the same four files with the
   history engine, re-equilibrate every increment and add whole-model finite
   differences of the ORIGINAL routine:

   ```bash
   resasm history --model Analysis.inp --odb Analysis.odb --material OTI_UMAT.obj \
       --request sensitivity_request.json --out results_checked --reequilibrate --verify fd
   ```

   Measured (1.6 s): `dU1/dE = -6.802721088435515e-09`,
   `dU1/dSIGY0 = -5.000000000000037e-04`, `dU1/dH = -1.2500000000000104e-05`,
   `dU1/dnu = -1.1e-17`: the closed form to about 1e-14. Its `run_report.txt`
   says `Tangent verified: yes: max relative error 1.46e-09 vs central FD of
   the ORIGINAL UMAT at 24 points` and `Derivative verified: yes: whole-model
   central FD of the ORIGINAL UMAT re-equilibrated in Python; worst
   nonzero-derivative error 2.66e-08 (plateau spread 1.67e-08)`.
3. **Replay against the ODB.** `run_report.txt` confirms that mesh, loads,
   supports, stresses, state and reactions of every increment were compared
   with the ODB; the GUI shows the stress difference, 4.1e-7 relative.
4. **Without Abaqus.** The export of this same ODB is committed under
   `tests/fixtures/presentation_j2/`. These offline tests replay it and assert
   the closed-form values (2 passed in 6.3 s on 2026-09-18):

   ```bash
   python -m pytest -q tests/integration/test_presentation_request.py -k "real_archived or compiled_offline"
   ```

## Run time

Measured on 2026-09-18: provider build 5.5 s; `resasm request` 1.1 s wall time,
including the ODB export through Abaqus Python; GUI Solve 1.0 s;
`resasm history --reequilibrate --verify fd` 1.6 s.

## Common problems

| Symptom | Cause and fix |
| --- | --- |
| `request failed: Category: odb_export` and exit 2 | Abaqus Python could not read the ODB: `abaqus` is not on `PATH`, or the ODB is incomplete. Pass the launcher with `--abaqus /path/to/abaqus`. `private/error_report.txt` holds the details. |
| `request failed: Category: output_directory` | The output folder already exists and is not empty. Results are never mixed; choose a new `--out`. |
| `request failed: Category: material_mapping` | No `Mapping.json` (or `OTI_UMAT.json`) beside the object, or it does not belong to this object (its SHA-256 or layouts differ). Copy the `.json` written by the same build, or pass `--mapping`. A valid mapping of another provider is not refused here: the request goes to the history engine, which names the reason on the first line (measured with the FCC provider: `outside the bounded presentation scope (mapping regular_source_hash does not match the pinned m3_j2 source: the material is umat_m6_fcc_oti, ...)`). |
| `--validate` stops with `Category: derivative_verification` | Measured on this ODB: its double-precision gate rejects the single-precision recorded displacements (`scaled error 1.731760e-08 >= 1.000000e-09` in `private/error_report.txt`). Use `resasm history --reequilibrate --verify fd` as shown above. |
| `umat-oti-provider build` fails | `gfortran` is missing or not on `PATH` ([INSTALL.md](../../docs/INSTALL.md)). |
| The Abaqus job ends with an error code although the `.sta` file says it completed | Judge the job by `THE ANALYSIS HAS COMPLETED SUCCESSFULLY` in the `.sta` file ([examples/cantilevers/README.md](../cantilevers/README.md) explains one such case). |

`resasm request` keeps models like this one on its bounded single-material
engine and hands every other readable model to the history engine, saying so
([Example 4](../replay_history/WALKTHROUGH.md), [Example 5](../cantilevers/WALKTHROUGH.md)).
All examples: [examples/README.md](../README.md).
