# GUI guide

Residual_Assembler has an optional graphical interface, a Streamlit web
application. It does not compute anything itself: every button builds the same
argument list you would type after `resasm`, calls the command line in the
same process, and shows the command, its exit code and its output verbatim. The
GUI therefore cannot show a number that the command line would not produce, and
every GUI run can be repeated in a terminal ([CLI guide](CLI_GUIDE.md)).

This guide walks through the screens in the order you are likely to use them:
first the **Sensitivity Request** screen, where most users spend their time,
then the assembly console (**Start here**, **1. Model** to **6. Backends**)
and **Advanced Replay**.

The behaviour described here was measured on 2026-09-18 by driving
`scripts/app.py` headless with Streamlit's own test harness
(`streamlit.testing`), and the server start by the health check of a real
server. The two screenshots are those in `docs/screenshots/`.

## Start the GUI

With the environment of [INSTALL.md](INSTALL.md) active (the `gui` extra
installs Streamlit), from the root of the Residual_Assembler checkout:

```bash
streamlit run scripts/app.py
```

Streamlit prints a local address (by default <http://localhost:8501>); open it
in a browser. Useful options:

```bash
streamlit run scripts/app.py --server.port 8502                  # another port
streamlit run scripts/app.py --server.headless true \
    --server.address 127.0.0.1 --server.port 8501                 # no browser pop-up
```

On a remote machine, forward the port with SSH (`ssh -L 8501:localhost:8501
<host>`) and open the local address. Stop the server with `Ctrl+C` in the
terminal.

Two things to set up before you work:

- **Where files are written.** The assembly tabs write into the sidebar's
  **Working directory**, and the request screen into its **Output directory**.
  Both default to `resasm_gui_workspace/` inside the checkout (ignored by git).
  For real work, point them to a folder outside the repository.
- **Where the tools are.** Reading an `.odb` needs Abaqus Python
  (`abaqus` on `PATH`, or its path under **Advanced**); the compiled material
  needs `gfortran`; the OTILib tabs need `PYOTI_PATH` and `OTILIB_ROOT` set in
  the shell that starts Streamlit ([INSTALL.md](INSTALL.md#4-install-otilib-for-the-direct-residual-paths)).

The page opens with the sidebar collapsed; open it with the arrow at the top
left.

## The screen layout

| Area | What it holds |
| --- | --- |
| Tabs along the top | **Sensitivity Request**, **Start here**, **1. Model**, **2. Requirements**, **3. Assemble**, **4. Sensitivity**, **5. Job**, **6. Backends**, **Advanced Replay** |
| Sidebar | **Where you are** (six steps, ticked only when the command really exited 0), the loaded model with **Clear model** or **Load the demo model**, the **Working directory** with **Create it**, the OTILib status with **What that means**, and a table of **Exit codes** |
| Each action | the `resasm` command it is about to run, the button, then a coloured result line: green `exit code 0 - success`, yellow `exit code 1 - ran, and the answer is negative`, red `exit code 2` (could not run) or `3` (OTILib missing), the command as run with its folder, and the full output |

A yellow result is an answer, not a crash: the output names the check that
failed or the ingredient that is missing.

## Sensitivity Request

![The Sensitivity Request screen after Solve](screenshots/resasm_solve.png)

This screen is `resasm request`: sensitivities of a finished analysis from four
inputs. Models inside the scope of the bounded single-material engine are
solved there; every other readable model is handed to the history engine, as on
the command line.

### 1. The four inputs

Each input has an **Upload** button (left) and a **path** field (right). Paths
are simpler for large files and keep the mapping discovery automatic.

| Field | What to give |
| --- | --- |
| **OTI_UMAT.obj** | the compiled material provider from `umat-oti-provider build`. Its mapping must sit beside it as `<object-stem>.json` or `Mapping.json`. If you upload the object instead of giving its path, upload the mapping too, under **Advanced** |
| **Analysis.inp** | the deck of the finished analysis |
| **Analysis.odb** | its output database |
| **sensitivity_request.json** | optional. If given, it is used exactly as written, and a caption says so. If empty, the controls below write the request for you |

### 2. Parameters, output and region

Under **3. Outputs and parameters** the screen reads the mapping and shows one
tick box per material parameter, all ticked (for the J2 material: `E`, `nu`,
`SIGY0`, `H`). Untick the ones you do not want. Then choose:

| Control | Choices |
| --- | --- |
| **output** | `displacement U1` to `U3`, `reaction force RF1` to `RF3`, `stress S11` to `S23`, `state variable SDV1` |
| **region** | `whole mesh`, or any node set (for U and RF) or element set (for S and SDV) of the `.inp`, read from the deck (for the one-element deck of Example 3: `node set LOADED`, `XZERO`, `YZERO`, `ZZERO`; `element set ALL`) |
| **summary over the region** | `mean`, `sum`, `max`, `L2`, `component` |
| **increments** | `last increment` or `every increment` |

These controls cover the core request format. For von Mises outputs,
`volume_mean`, several outputs at once, weighted shares or a full-field
archive, write a request file ([format](CLI_GUIDE.md#the-request-file)) and
give it in the fourth field.

### 3. Output directory and Advanced

**Output directory** must be a folder that does not exist yet, or is empty.
Under **Advanced**:

| Field | Meaning | CLI |
| --- | --- | --- |
| **Mapping.json path** / **Mapping.json** upload | the mapping, when it is not beside the object | `--mapping` |
| **Abaqus executable** | default `abaqus`; the launcher used to read the ODB | `--abaqus` |
| **Independent finite-difference validation** | off by default | `--validate` |

### 4. Solve

**Solve** is enabled once the object, the deck and the ODB are given and there
is a request (a file, or at least one ticked parameter). It runs:

```bash
resasm request --out <Output directory> --abaqus <Abaqus executable> \
    --material <object> --model <deck> --odb <ODB> --request <request file> [--mapping ...] [--validate]
```

When no request file is given, the choices are first written to a
`sensitivity_request.json` in a temporary folder. The request actually used is
recorded under `request` in `sensitivity_results.json`, so any GUI run can be
repeated on the command line.

### 5. What appears after Solve

- A green line, for example `Executed: 1 scalar results. Independent
  validation: not run.`
- Three download buttons: `sensitivity_results.json`,
  `sensitivity_tables.csv` and `run_report.txt`, the public results. Everything
  else stays in the output directory, under `private/`.
- For a run of the bounded engine, the **full field** of the chosen region:
  one row per node (U) or per integration point (S, SDV) with the value, its
  derivative with respect to every ticked parameter and the **governing
  parameter** (the largest `|p dy/dp|` at that location, with `p` the `.inp`
  constant at the mapping's PROPS index); and two figures, **stress reproduced
  vs. the .odb** (largest difference between replayed and recorded stress,
  relative) and **solution sensitivity ‖du/dp‖**. The table's own toolbar
  downloads it as CSV. Reaction derivatives are not stored per node, so a
  reaction output shows only the reduced value.
- For a run handed to the history engine, the downloads only. The full field of
  such a run is the `fields.npz` in the output directory, written when the
  request file asks for `"full_field": true`
  ([Example 5](../examples/cantilevers/WALKTHROUGH.md) shows how to read it).

### Measured runs

On the one-element J2 analysis of
[Example 3](../examples/presentation_request/WALKTHROUGH.md), all four
parameters ticked, **output** `displacement U1`, **region** `node set LOADED`,
**summary** `mean`, **increments** `last increment`: Solve took 1.0 s, reported
`Executed: 1 scalar results. Independent validation: not run.`, and the table
had four rows (nodes 2, 3, 6, 7), each with value `0.026429`,
`d/dE = -6.802720e-09`, `d/dnu = 1.172950e-09`, `d/dSIGY0 = -5.000000e-04`,
`d/dH = -1.250000e-05` and governing parameter `SIGY0`; **stress reproduced
vs. the .odb** `4.1e-07`, **‖du/dp‖** `0.00422`. This is the run in the
screenshot above. The numbers equal those of `resasm request` on the same
files, and the uniaxial closed form within the precision of the ODB.

The same screen with **output** `state variable SDV1` and **region** `whole
mesh` gave eight rows, one per integration point, each with `EQPLAS = 0.025`,
`d/dSIGY0 = -5.000000e-04` and `d/dH = -1.250000e-05`:

![Equivalent plastic strain and its derivatives at every integration point](screenshots/resasm_solve_sdv_all_points.png)

On the full-size J2 cantilever of
[Example 5](../examples/cantilevers/WALKTHROUGH.md), with
`examples/cantilevers/j2_request.json` given as the request file, Solve took
24.5 s and reported `Executed: 240 scalar results. Independent validation:
not run.` with the three downloads.

### An output directory can be used once

A second Solve into the same **Output directory** is refused:

```text
request failed: Category: output_directory
Action: use a new or empty writable output directory; refusing stale public results.
Private diagnostics: unavailable (output directory could not be prepared or written).
```

This is intended. The request code refuses any existing, non-empty output
folder (`run_request` in `residual_core/replay/presentation.py`; the history
engine has the same guard), and the offline test suite checks that it does.
Without it, a `sensitivity_results.json` from an earlier run could sit next to
files of a later one, and results from two different requests could be taken
for one. The default `resasm_gui_workspace/request` therefore works for the
first Solve only. For each new Solve, type a new folder name (for example
`.../request_U1_loaded`, `.../request_sdv1`), or delete the old folder
yourself if you no longer need it. Nothing is ever overwritten for you.

### Other refusals

| Message (category) | What to do |
| --- | --- |
| `odb_export` | Abaqus Python could not read the ODB. Check **Abaqus executable** under **Advanced** and that the ODB is complete. |
| `material_mapping` | No mapping beside the object, or it belongs to another object. Keep the `.json` from the same build beside the object, or give it under **Advanced**. |
| `request_schema` | The request file does not match the format ([CLI guide](CLI_GUIDE.md#the-request-file)). |
| `derivative_verification` | The optional validation disagreed or refused single-precision data; see [Example 3](../examples/presentation_request/WALKTHROUGH.md#common-problems). |
| `Mapping.json could not be read` above the tick boxes | The file next to the object is not a completed mapping. |

Each refusal is also written to `run_report.txt` in the output folder, with
details in `private/error_report.txt`.

## Start here

A one-button tour that needs no Abaqus and no OTILib. **Run the demo** loads
the unit-cube C3D8 of
[Example 2](../residual_core/examples/minimal_c3d8_stress_driven/WALKTHROUGH.md)
(with its stress field) and runs four commands, showing each with its output:

| Step | Command | Measured result |
| --- | --- | --- |
| 1 | `resasm inspect <model>` | exit 0 |
| 2 | `resasm requirements <model> --mode stress-driven --fields <fields>` | exit 0, `Ready to assemble in stress-driven mode.` |
| 3 | `resasm assemble ... --out demo_R.npy` | exit 0, `ndof=24  ||R||=7.071068e+01  max|R|=2.500000e+01` |
| 4 | `resasm verify <model> --fields <fields>` | exit 1: the unloaded cube is not in equilibrium, as explained in Example 2 |

A **Download demo_R.npy** button (320 bytes) follows. After the demo the model
stays loaded, so every other tab is unlocked. The tab also explains the exit
codes.

## 1. Model

Everything in tabs 2 to 4 needs a model. Load one here:

- **Shipped examples**: the **Example model** list offers the neutral-format
  models under `residual_core/examples/` (`minimal_beam`,
  `minimal_c3d8_stress_driven`, `minimal_mixed`,
  `minimal_nonlinear_spring_sensitivity`, `minimal_truss`). Choosing one loads
  it at once; its README summary and files are shown, and a `fields.json` beside
  it is filled into the `--fields` boxes of tabs 2 and 3.
- **Your own model**: type an absolute path to an `.inp` or `.json` and press
  **Use this path**, or upload a file and press **Use uploaded file** (it is
  saved into the working directory; an `.inp` that includes other files needs
  them too).

Then **Run inspect** (tick `--detail` for the per-element backend selection)
and **List modes**.

CLI: `resasm inspect [--detail] <model>`, `resasm modes`.

## 2. Requirements

- **requirements --mode**: choose the **Mode**, optionally `--fields` and
  `--subroutine`, and press **Run requirements**. It names the single missing
  ingredient, or says `Ready to assemble in <mode> mode.`
- **doctor**: press **Run doctor** for inspection plus a readiness line per
  mode. Tick `--write-config-template` to also write a config template to the
  path shown; a download button appears.

CLI: `resasm requirements <model> --mode <mode> [--fields ...] [--subroutine ...]`,
`resasm doctor <model> [--write-config-template <file>]`.

## 3. Assemble

- **assemble**: choose `--mode`, fill `--fields` or `--subroutine` as the mode
  needs, tick `--tangent` to assemble the tangent too, and set `--out`
  (default `R.npy`, in the working directory). **Run assemble**, then
  **Download R.npy**. On the stress-driven cube the measured output was
  `assembled mode 'stress-driven': ndof=24  ||R||=7.071068e+01  max|R|=2.500000e+01`.
  If the mode is not runnable, the result is red (exit 2) and shows the
  requirements report instead of a residual.
- **verify**: **Run verify** compares the free residual with the tolerance and
  lists the reactions (exit 1 on the unloaded demo cube, as expected).

The tab does not pass the global `--config` option. For a model that needs a
config (the finite-strain example), assemble on the command line
([Example 7](../examples/finite_strain_c3d8/WALKTHROUGH.md)).

CLI: `resasm assemble <model> --mode <mode> [--fields ...] [--subroutine ...] [--tangent] [--out R.npy]`,
`resasm verify <model> [--fields ...]`.

## 4. Sensitivity

The tab first states whether OTILib is available. Then:

- `--params` is filled with the `params.json` beside the model, if there is
  one; or list parameters one per line under `--param`.
- Choose `--mode`, `--order` (0 leaves it to the params file), `--backend`,
  tick `--no-fd` to skip the finite-difference cross-check, and give `--out` to
  save the package (three download buttons appear).
- Press **Run sensitivity**.

Measured on `minimal_nonlinear_spring_sensitivity` (exit 0):
`d^1/e1 = -3.333333e-01   [analytic -3.333333e-01, rel 0.00e+00]   [FD -3.333330e-01, rel 1.00e-06]`.
On the finite-strain model of Example 7, loaded by path, the same button gave
`+1.872427e-01` and `+3.398245e-02` with finite-difference agreement 2.49e-10
and 6.79e-10. The finite-difference column is computed by a different method
(re-solving the model at perturbed parameters). If any derivative disagrees
with it by `1e-4` or more, the command says `finite-difference check FAILED`
and the result is yellow (exit 1). A backend that cannot produce a derivative
refuses (red, exit 2) instead of printing zeros: `--backend dual1` on the
finite-strain model names the reason.

CLI: `resasm sensitivity <model> [--params ...] [--param ...] [--mode ...] [--order N] [--backend ...] [--no-fd] [--out PREFIX]`.

## 5. Job

The `resasm.yml` workflow in five steps, all relative to the working
directory. One **Job folder** name drives the tab.

| Step | Button | CLI |
| --- | --- | --- |
| a. Start from a template | choose **--template** (`blackbox`, `blackbox-order2`, `cpp`, `fortran`, `python`), optionally **--force**, press **Copy template**; the current `resasm.yml` is shown below it | `resasm init --template <t> --out <job>` |
| b. Or write a recipe from a model | fill `--solution`, `--material`, `--param`, `--name`, `--out`; **Run inspect-model**, **Write recipe** (give the model as an absolute path) | `resasm inspect-model ...`, `resasm init-assembly ...` |
| c. Check the config | **Run check** | `resasm check <job>/resasm.yml` |
| d. Run the job | **Run job** | `resasm run <job>/resasm.yml` |
| e. Read the report | **Read report**; `public/summary.md` opens below with a download button | `resasm report <job>/resasm.yml` |

Measured with the `python` template (OTILib available): **Copy template**,
**Run check**, **Run job** and **Read report** all exited 0, and
`<job>/resasm_output/public/summary.md` opened below step e.

Step e gives the job's `resasm.yml` to `resasm report`, which reads the
folder that job writes to: `<job>/resasm_output`, or the folder named by
`output: dir:` in `resasm.yml`. Measured with the `blackbox-order2` template
and `output: dir: results_here` added before step d: **Read report** exited 0
and opened `<job>/results_here/public/summary.md`.

The interactive `resasm init` wizard (without `--template`) reads answers from
a terminal and is available on the command line only.

## 6. Backends

- **Audit all backends** runs `resasm backends`: every registered element
  backend with its element types, modes, tangent, status and limitations.
- **template**: choose `--formulation` or `--material`, pick a **Name**, and
  press **Print contract** (`resasm template --formulation <name>`).

## Advanced Replay

The bounded J2 replay of a record or a small model with total-history
sensitivities (`resasm replay`), used by
[Example 6](../examples/bounded_j2_c3d8/WALKTHROUGH.md).

| Field | Default | CLI |
| --- | --- | --- |
| **Record or model JSON** | `examples/bounded_j2_c3d8/model.json` | `RECORD` |
| **Compiled provider object** | empty | `--object` |
| **Provider contract JSON** | empty | `--contract` |
| **Replay output directory** | `resasm_gui_workspace/replay` | `--out` |
| **Solve model** | ticked | `--solve` |
| **Verify with ORIGINAL finite differences** | ticked | `--verify` |

Fill the object and the contract (the `.obj` and `.json` written by
`umat-oti-provider build`) and press **Run replay**. Measured with a freshly
built J2 provider: exit 0 in 1.60 s,
`replay: 7 increments, 8 IPs; total-history du/dp; verified=True`.

## Which tab is which command

| Tab | Commands |
| --- | --- |
| Sensitivity Request | `resasm request` (which may hand the model to `resasm history`) |
| Start here | `inspect`, `requirements`, `assemble`, `verify` |
| 1. Model | `inspect`, `modes` |
| 2. Requirements | `requirements`, `doctor` |
| 3. Assemble | `assemble`, `verify` |
| 4. Sensitivity | `sensitivity` |
| 5. Job | `init`, `inspect-model`, `init-assembly`, `check`, `run`, `report` |
| 6. Backends | `backends`, `template` |
| Advanced Replay | `replay` |

Not in the GUI: `resasm history` with `--fields`, `--reequilibrate` or
`--verify` (use the command line; [Example 4](../examples/replay_history/WALKTHROUGH.md)),
the global `--config` option, and the interactive `init` wizard.

## Troubleshooting

| Symptom | Cause and fix |
| --- | --- |
| `streamlit: command not found` | Install the `gui` extra (`pip install -e "./Residual_Assembler[gui]"`) and activate the environment. |
| **Solve** stays grey | The object, the deck and the ODB must all be given, and there must be a request: a file, or at least one ticked parameter. |
| No tick boxes appear | The mapping was not found: keep it beside the object, or give it under **Advanced**. |
| `request failed: Category: output_directory` | The output folder is not new or empty; see [above](#an-output-directory-can-be-used-once). |
| A tab says **Locked** | Load a model on **1. Model**, or press **Load the demo model**. |
| The sidebar says `OTILib: not installed` | Set `PYOTI_PATH` and `OTILIB_ROOT` in the shell, then restart Streamlit. |
| The **Where you are** list lags by one step | The sidebar is drawn before the tab that ran the command; it catches up at the next interaction. |

Further reference on the request screen and its tests:
[GUI.md](GUI.md). The worked examples: [examples/README.md](../examples/README.md).
