# Claude agent G: the slide-18/42 Solve screen

Branch `claude/G-2026-09-18`, from the snapshot `d8ea804` of Copilot's
`imq-ra-recovery` working tree. This file lists **every edit to a file Copilot
also edits**, for the merge. New files follow.

## Edits to files Copilot owns or edits

### `residual_core/app/streamlit_app.py` (`_tab_request` only)

The screen already took the object, the .inp, the .odb and a hand-written
`sensitivity_request.json`, and it ran `resasm request`. The slide also shows
parameter ticks, an output, a region, a Solve button and the full-field result.
Four hunks, 13 lines added and 2 changed:

1. `from residual_core.app import request_screen` at the top of `_tab_request`.
2. After the four input rows: `generated = request_screen.render_outputs_and_parameters(...)`
   draws section "3. Outputs and parameters" (ticks from Mapping.json, output,
   region, summary, increments) and returns the request those choices stand for.
3. The button: label `"Run sensitivity request"` became `"Solve"` (key
   `btn_request_run` unchanged, `type="primary"` added). `disabled=` no longer
   requires a request file when the choices make a request. Inside the
   temporary-directory block, when no request file was given, the generated
   request is written with `request_screen.write_request` and passed as
   `--request`. Everything else in the run path is unchanged.
4. After the three downloads: `request_screen.render_full_field(st, completed)`.

Copilot's AppTest (`tests/integration/test_presentation_request.py::test_gui_shared_service_with_compiled_core`)
passes unchanged. It still sees exactly three download buttons: the
full-field table uses the dataframe's own CSV download instead of a fourth
button, for this reason.

### `scripts/check_presentation_browser.py`

One line: it clicks `"Solve"` instead of `"Run sensitivity request"`. Run on
2026-09-18 against `../imq_abaqus/recovery_presentation/imqrp_reference/collaborator`
with output under `../imq_abaqus/claude_G/`: `passed: true`, three downloads,
no horizontal overflow at 390 px, and the server stopped.

### `docs/PRESENTATION_INTERFACE.md`

Four sentences after the GUI paragraph: the slide-18 widgets, the Solve
button, the full-field table, and a link to docs/GUI.md.

### `pyproject.toml`

One marker line: `gui` (browser workflow tests, deselected unless `-m gui`).

### `tests/framework/test_gui_is_a_thin_cli_front_end.py` (pre-existing, not Copilot's)

`test_the_app_does_not_import_the_assembler_directly` now walks
`residual_core/app/request_screen.py` too, with the same rule. The allow-list
gained that module and `residual_core.io.abaqus_inp_parser`, which is only
used to list the .inp's sets for the region dropdown. In the first draft,
`request_screen` called `presentation.scalar_results` to show per-node
reactions. The test caught it, so that call was removed: a reaction output now
shows only the reduced value the request computed.

## New files

- `residual_core/app/request_screen.py`: the tick boxes (from Mapping.json),
  the output and region choices, the request they write, and the full-field
  view. The view reads `private/result.json` and `private/fields.json` and
  imports nothing from the engine.
- `tests/gui/test_imqcam_solve_screen.py` (AppTest, offline suite, 6 tests),
  `tests/gui/test_imqcam_solve_screen_browser.py` (`-m gui`, abaqus, 1 test),
  `tests/gui/gui_helpers.py`, `tests/gui/conftest.py`.
- `docs/GUI.md`, `docs/screenshots/resasm_solve.png`,
  `docs/screenshots/resasm_solve_sdv_all_points.png`.

## Commands and results (2026-09-18)

Environment: `PYTHONPATH=<this worktree>:<claude-umat-G>/src:<otilib build_py311>`,
`UMAT_OTI_REPO=<claude-umat-G>`, Python 3.11 venv, gfortran 9.4, Abaqus
2021.HF5.

```text
python -m pytest -q -p no:cacheprovider tests/gui/test_imqcam_solve_screen.py \
  tests/integration/test_presentation_request.py tests/framework/test_gui_is_a_thin_cli_front_end.py
    -> 45 passed
python -m pytest -q -p no:cacheprovider -m gui tests/gui
    -> 1 passed, 6 deselected (real ODB through abaqus python, two Solves, 41 s)
```

Real-ODB run (inputs: provider built by `python -m umat_oti.provider.collaborator`
for m3_j2; `Analysis.inp` SHA-256 1cde25bd...; `Analysis.odb` 54745b97...).
The choices were E, nu, SIGY0 and H; U1; node set LOADED; mean. Value
0.026428570970892906. dE -6.802719961317867e-09 (reference -300/E² =
-6.8027e-09), dSIGY0 -4.999999999999996e-04 (-1/H), dH -1.249999988950807e-05
(-(300-SIGY0)/H²), dnu 1.1729503437354757e-09 (reference 0; see docs/GUI.md).
`resasm request` with the recorded request gave identical results. Stress was
reproduced to 4.1e-07.

`resasm request --validate` on the real ODB is refused by the engine's own
double-precision gate ("ORIGINAL primal displacement: scaled error 1.73e-08 >=
1e-09"), as docs/PRESENTATION_INTERFACE.md says it may be. Independent
validation of the derivatives on a float32 ODB is therefore not established by
this screen.

## Not done here

- The slide's output "von Mises stress σ_vM" and its crystal-plasticity
  cantilever: `resasm request` supports neither. The screen offers exactly
  what the request engine supports.
