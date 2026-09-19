# The Residual_Assembler GUI

One Streamlit application over the `resasm` command line. It opens on
**Sensitivity Request**, the screen that computes the sensitivities of a
finished analysis. The other tabs (**Start here**, **1. Model** to
**6. Backends**, **Advanced Replay**) are the assembly console. This page is
reference material on the Sensitivity Request screen and its tests; the
step-by-step guide to every screen is [GUI_GUIDE.md](GUI_GUIDE.md).

Every button runs the real command in-process through
`residual_core.ui.cli.main` and shows its exit code, so the GUI cannot report
a number the CLI would refuse to produce
(`tests/framework/test_gui_is_a_thin_cli_front_end.py` pins this; it also
covers the Solve screen's helper module).

## Launch

```bash
pip install -e ".[gui]"
streamlit run scripts/app.py
```

Reading an `.odb` needs licensed Abaqus Python on `PATH` (`abaqus python`).
Running the compiled material needs the matching Fortran runtime (gfortran on
the verified Linux machine). The material source is never needed.

## Sensitivity Request: point, tick, choose, Solve

![The Solve screen after a run on the saved J2 analysis](screenshots/resasm_solve.png)

You point at the compiled OTI object from the provider build, and at the
saved analysis: its `.odb` and the `.inp` that produced it. You tick the
parameters and choose the output and the region, then press **Solve**. The
residual equation is assembled and solved over the recorded history, and the
screen shows the sensitivities in the chosen region.

| Field | What it is |
| --- | --- |
| `OTI_UMAT.obj` | the compiled OTI object from the provider build (upload or path). Its `Mapping.json` must sit beside it, or be chosen under Advanced |
| `Analysis.inp`, `Analysis.odb` | the saved analysis and the input that produced it |
| `sensitivity_request.json` | optional. When it is given, it is used as written. When it is not, the choices below write it |
| parameter tick boxes | one per parameter named in `Mapping.json`, all ticked at first |
| output | a displacement component (U1-U3), a reaction component (RF1-RF3), a stress component (S11-S23) or state variable SDV1: exactly what `resasm request` supports |
| region | the whole mesh, or any node set (for U, RF) or element set (for S, SDV) of the `.inp` |
| summary over the region | mean, sum, max, L2 or component, the reductions of `resasm request` |
| increments | the last increment or every increment |

**What it calls.** `resasm request` (`residual_core/ui/cmd_request.py`, whose
`route_request` runs `run_request` in `residual_core/replay/presentation.py`
for models inside the bounded engine's scope and the history engine for the
others). The parameter ticks
and the output and region choices only write the `sensitivity_request.json`
that a user would otherwise write by hand
(`residual_core/app/request_screen.py`). The run is exactly:

```bash
resasm request --model Analysis.inp --odb Analysis.odb --material OTI_UMAT.obj \
  --request sensitivity_request.json --out results
```

The request actually used is recorded in `results/sensitivity_results.json`
(`request`), so the command line can repeat any GUI run.

**What it writes.** The output directory (it must be new or empty) holds the
three public files, which the screen offers for download:
`sensitivity_results.json`, `sensitivity_tables.csv` and `run_report.txt`.
It also holds `private/`. For the bounded engine that is the full history
arrays in `result.json`, the strict ODB export in `fields.json`, the export
command and log, and the link shim; for the history engine, the files listed
in [REPLAY_HISTORY.md](REPLAY_HISTORY.md#outputs-and-report).

**What it shows after Solve.**

- The full field in the chosen region: one row per node (U) or per
  integration point (S, SDV), with the value and its derivative with respect
  to every ticked parameter. These are the arrays the request stored in
  `private/result.json` (`u`/`du_dp`, `stress`/`dsigma_dp`,
  `state`/`dstate_dp`). The table's own toolbar downloads it as CSV.
  Reaction derivatives are not stored per node, so a reaction output shows
  only the reduced value.
- "governing parameter": the ticked parameter with the largest |p·∂y/∂p| at
  that location. p is the `.inp`'s USER MATERIAL constant at the
  `Mapping.json` PROPS index.
- "stress reproduced vs. the .odb": max |replayed stress − ODB stress| /
  max(|ODB stress|, 1) at the last output increment.
- "solution sensitivity ‖du/dp‖": the Frobenius norm of du/dp over every
  degree of freedom and the ticked parameters.

**Measured (2026-09-18), real ODB.** The object and `Mapping.json` came from
the provider build of UMAT-OTI (`umat_oti.provider.collaborator`, the function
behind its GUI Build button) for the pinned m3_j2 source. The saved analysis
is the one-element example `examples/presentation_request/Analysis.inp` (one
C3D8, four increments, E=210000, nu=0.3, SIGY0=250, H=2000) and its Abaqus
2021.HF5 ODB (SHA-256 54745b97...). The choices were all four parameters,
displacement U1, node set LOADED, mean. Solve reported 1 scalar result with
value 0.026428570970892906 and derivatives

| parameter | GUI (= `resasm request`) | uniaxial reference |
| --- | --- | --- |
| E | -6.802719961317867e-09 | -300/E² = -6.8027e-09 |
| nu | 1.1729503437354757e-09 | 0 (U1 of a uniaxial stress state does not depend on nu). A double-precision forward solve of the same model through the same object gives \|dU1/dnu\| ≤ 1.3e-16. The replay starts from the recorded displacements, which are not re-equilibrated (scaled free residual 4.1e-07), and gives 1.17e-09. That is 1.3e-8 of U1 after scaling by nu. |
| SIGY0 | -4.999999999999996e-04 | -1/H = -5.0e-04 |
| H | -1.249999988950807e-05 | -(300-SIGY0)/H² = -1.25e-05 |

The same request through `resasm request` gave identical numbers. Stress was
reproduced to 4.1e-07, ‖du/dp‖ was 0.00422, and the maximum scaled free
residual was 4.12e-07.

A second Solve on the same analysis chose state variable SDV1 (equivalent
plastic strain) over the whole mesh. The table has one row per integration
point, all eight of them. Each has EQPLAS 2.500000e-2 = (300-SIGY0)/H,
dEQPLAS/dSIGY0 = -5.000000e-4 = -1/H and dEQPLAS/dH = -1.250000e-5 =
-(300-SIGY0)/H². dE is at most 9.3e-23 and dnu at most 3.1e-17, both zero to
rounding as the uniaxial reference says:

![Equivalent plastic strain and its derivatives at every integration point](screenshots/resasm_solve_sdv_all_points.png)

## Tests

| Test | What it drives | Run |
| --- | --- | --- |
| `tests/gui/test_solve_screen.py` | the screen through `streamlit.testing` (AppTest), with a freshly built OTI object and its Mapping.json. The ODB export is replaced by the genuine export of the same ODB (`tests/fixtures/presentation_j2/fields.json`). Checks the ticks, output, region and Solve for U on node set LOADED, S11 and SDV1 at every integration point, and RF1 on XZERO against the uniaxial references and against `resasm request`. Also checks that a supplied request file still wins and that there is nothing to tick without a mapping | the offline suite |
| `tests/gui/test_solve_screen_browser.py` | the provider hand-over (object and mapping from the build) consumed by the Sensitivity Request screen in headless Chromium, on the real ODB (licensed Abaqus Python): U1 on node set LOADED, then SDV1 on the whole mesh. Writes `docs/screenshots/resasm_solve.png` and `docs/screenshots/resasm_solve_sdv_all_points.png` | `python -m pytest -m gui tests/gui` (also marked `abaqus`) |
| `tests/integration/test_presentation_request.py` | the request-file route of `resasm request` | the offline suite |

`scripts/check_presentation_browser.py` drives the request-file route in a
browser, at desktop and phone widths, including the **Solve** button and the
three downloads.

## Limits

- The screen's own controls write the core request format: one output among
  U1 to U3, RF1 to RF3, S11 to S23 and SDV1, a region, a reduction among
  mean, sum, max, L2 and component, and the last or every increment. For von
  Mises outputs, volume-weighted means, several outputs, weighted shares or a
  full field, give a request file ([CLI_GUIDE.md](CLI_GUIDE.md#the-request-file)).
- Solve runs `resasm request`, so a model outside the bounded single-material
  engine is handed to the history engine exactly as on the command line
  ([REPLAY_HISTORY.md](REPLAY_HISTORY.md#how-resasm-request-chooses-this-engine)).
  Measured on 2026-09-18: the full-size J2 cantilever of
  [examples/cantilevers](../examples/cantilevers/WALKTHROUGH.md), with its
  request file, solved from this screen in 24.5 s with 240 scalar results.
  What the screen shows after such a run is described in
  [GUI_GUIDE.md](GUI_GUIDE.md#5-what-appears-after-solve).
- Each Solve needs a new or empty output directory; a used one is refused so
  that results of two runs are never mixed
  ([GUI_GUIDE.md](GUI_GUIDE.md#an-output-directory-can-be-used-once)).
