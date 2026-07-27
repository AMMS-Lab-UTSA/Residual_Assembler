# Definition deck — slide-by-slide outline

Twelve slides defining Residual_Assembler: a privacy-preserving residual-sensitivity
toolkit. This is an outline, not a deck. Every claim below is grounded in the
current code; slide 11 states plainly what is verified and what is pending.

---

## Slide 1 — What problem are we solving?

- A user has a working model. They want parameter sensitivities `du/da` from it.
- Their code, mesh, and material models are private. They cannot ship them out, and
  they cannot let an outside tool read them.
- Finite differences need one extra nonlinear solve per parameter, and the step size
  is a guess.
- We only ever need to *call* the residual, never to read it. That is enough to get
  exact, arbitrary-order sensitivities.
- One residual function. One small `resasm.yml`. One command. A sensitivity package out.

**visual:** box diagram — user's private residual (locked box) on the left, `resasm`
in the middle, shareable public report on the right; a dashed "privacy boundary" line
between left and middle that no source/mesh/state arrow crosses.

---

## Slide 2 — The residual sensitivity equation

- At the converged real solution, `R(u, a) = 0`. Differentiate w.r.t. the parameters
  and you get the linear system we solve, once per order `p`:

  ```
  T U^(p) = -R^(p)
  ```

- `T` = tangent / Jacobian `dR/du` at the converged real solution. The user supplies it.
- `R^(p)` = p-th order residual coefficient matrix, read off the OTI evaluation of
  the user's residual.
- `U^(p)` = p-th order solution sensitivity coefficients. This is the answer.
- Same `T` at every order — factor once, reuse. Only the right-hand side changes.
- Solved on the free DOFs only; prescribed DOFs get zero sensitivity
  (`resasm_user/oti_global.py`).

**visual:** the equation large in the centre, with each of the three symbols
labelled by who produces it: `T` = user, `R^(p)` = OTI evaluation, `U^(p)` = output.

---

## Slide 3 — What the user provides

- A **converged solution vector** `u` (`solution.file` as `.npy`, or `solution.values`).
- **Parameter names and real values** — `parameters: {k: 2.0, f: 16.0}`. These are the
  design variables that get seeded.
- A **residual evaluator** — Python function, executable, or element module.
- A **tangent** `T = dR/du` at `u` — a Python `tangent()`, a saved `tangent.npz`, or
  returned in the black-box response.
- A **derivative order** — `sensitivity.order` (≥ 1).
- **No mesh, no element data, no material model** is needed when the residual returns
  the *global* residual (`docs/user_input_contract.md`).

**visual:** a five-row "inputs" box (u, parameters, residual, tangent, order) beside a
crossed-out "not required" box (mesh, elements, materials, source code, state).

---

## Slide 4 — The residual provider contract

- Three provider paths, one contract, identical output package
  (`docs/which_path_should_i_use.md`).
- **Python** (`residual.type: python`): implement
  `residual(u, params, state=None, time=None) -> R` using ordinary arithmetic. We call
  it twice — once with floats, once with OTI scalars. Template:
  `templates/user_python_residual/`.
- **Black-box executable** (`residual.type: executable`): your program reads
  `request.json` (`u`, `parameters`, `seed_directions`, `basis_count`,
  `truncation_order`, `order`, `direction_map`, `u_star_coefficients`, `state`, `time`)
  and writes `response.npz` with arrays `R_order_<p>` of shape `(unknowns, N^(p))`,
  plus an optional `tangent`. Template: `templates/user_blackbox_residual/`.
- **Compiled C++/Fortran**: the same request/response contract, wrapped around
  existing compiled code. Templates: `templates/user_cpp_residual/`,
  `templates/user_fortran_residual/`; headers in `partner_kit/`.
- Only the Python path needs OTILib on our side. The executable paths produce `R^(p)`
  themselves, so we never link to the user's code at all.

**visual:** three lanes into one funnel — "Python callable", "executable +
request/response JSON/NPZ", "compiled C++/Fortran wrapper" → one box labelled
`R^(p)`, `T`.

---

## Slide 5 — The OTI execution loop

- OTI/HYPAD in one screen (`resasm_user/oti_global.py::solve_python`):
  1. **Define the number system.** `m` = number of parameters, truncation order `q`.
     One OTI number carries every parameter direction at once.
  2. **Seed.** `a_i* = a_i + e_i` for every parameter — all `m` in one shot, not one
     run per parameter.
  3. **Evaluate.** `R* = residual(u*, a*)` — the *user's* function, called with OTI
     scalars.
  4. **Extract.** `R^(p)` = the order-`p` coefficients of `R*`, one column per
     direction (including mixed directions like `e1*e2`).
  5. **Solve and re-inject.** `U^(p) = -T_ff^-1 R^(p)` on the free DOFs, then write
     `U^(p)` back into `u*` **before** order `p+1`.
- Loop `p = 1 .. q`. One residual evaluation per order — not one per parameter.
- The order-`p` direction list and each direction's `recovery_factor` (which converts
  an OTI coefficient into the true partial derivative) come from
  `residual_core/algebra/otilib_adapter.py::order_directions`.

**visual:** a 5-box loop (seed → evaluate → extract `R^(p)` → solve `T U^(p) = -R^(p)`
→ inject into `u*`) with the arrow from box 5 back into box 2 labelled "p → p+1".

---

## Slide 6 — The output package

- One directory per run (default `resasm_output/`), split `private/` and `public/`
  (`resasm_user/output.py`, `docs/user_output_contract.md`).
- **private/** — `metadata.json`, `parameter_map.json`, `dof_map.json`,
  `residual_real.npz`, `tangent.npz`, `rhs_order<p>.npz` (`R^(p)` and `-R^(p)`),
  `solution_sensitivities_order<p>.npz` (`U^(p)`), `validation_full.json`.
- **public/** — `summary.md`, `parameter_ranking.csv` (ranked by `||du/da||`),
  `sensitivity_norms.csv` (per-direction solution and RHS norms),
  `validation_summary.json`, `timing.json`.
- `resasm report resasm_output/` prints the headline: residual norm, tangent source,
  parameters, orders solved, validation status.
- Machine-readable: `from resasm_user import read_report`.

**visual:** a two-column file tree, `private/` shaded red ("stays local"), `public/`
shaded green ("safe to send"), with an arrow labelled "share" leaving only the green
column.

---

## Slide 7 — Privacy model

- The tool calls the residual; it never reads it. On the black-box path it only
  exchanges two files with an executable it knows nothing about.
- **What is never required:** source code, mesh, element connectivity, material
  models, state variables.
- **What `public/` never contains:** source, mesh, the full residual vector, the full
  tangent matrix, state variables, or private parameter *values* — only parameter
  **names** and sensitivity **norms**. Enforced where the outputs are written.
- **What `private/` holds:** every full array. It stays on the user's disk unless they
  choose to send it.
- Honesty rule in the report: on the black-box path the real residual is never
  exposed, so `residual_free_norm` is written as `null` and `summary.md` says
  "Equilibrium was NOT verified by this run" — never a fabricated `0.0`.

**visual:** the privacy boundary from slide 1, redrawn as a table — rows = artefacts
(source, mesh, residual vector, tangent, state, parameter values, parameter names,
norms), columns = "private / public", with only the last two ticked under public.

---

## Slide 8 — Verification ladder

- `resasm check resasm.yml` — a step-by-step readiness report. Each line is `[ok]`,
  `[warn]` or `[fail]`; it stops at the first actionable blocker and shows the exact
  YAML to add. It parses the config, loads the solution, evaluates the residual,
  loads the tangent, confirms OTILib, checks `||R_free||` (warns if `u` may not be
  converged), then does a real order-1 RHS generation and solve.
- **RHS finite-difference cross-check** (`validation.rhs_finite_difference_check`):
  central-differences `d(residual)/d(a_i)` at the fixed `u`, solves with the **same**
  tangent, compares to `U^(1)`. Pass threshold: max relative error < 1e-4.
- **What that check does NOT prove** — and we say so in the report: it does not
  re-solve the nonlinear problem, so a wrong tangent cancels out on both sides. It is
  not a solution-level FD validation. `validation.solution_finite_difference_solver`
  is a reserved config field and is **not implemented**; setting it prints exactly that.
- The FE framework side has its own generic ladder, Levels 0–7 (zero-field,
  rigid-body, patch test, FD tangent, solver comparison, material replay, full
  residual replay) — `residual_core/core/verification.py`.

**visual:** a ladder graphic, rungs 0–7, with the rungs reachable offline shaded solid
and the Abaqus-dependent rungs (5, 6, 7) drawn as dashed outlines.

---

## Slide 9 — Example: one-line spring residual

- `templates/user_python_residual/` — the whole model is two functions:

  ```python
  def residual(u, params, state=None, time=None):
      return [params["k"] * u[0] ** 3 - params["f"]]     # R = k u^3 - f

  def tangent(u, params, state=None, time=None):
      return [[3.0 * params["k"] * u[0] ** 2]]
  ```

- Config: `parameters: {k: 2.0, f: 16.0}`, `solution.npy` = `[2.0]`,
  `sensitivity: {order: 2, backend: otilib}`.
- Analytic answers: `du/dk = -u/(3k) = -1/3`, `du/df = 1/(3k u^2) = 1/24`. The run
  reproduces them, and the RHS FD cross-check confirms order 1.
- Order 2 comes from the same single run — the OTI number already carries `e1^2`,
  `e1*e2`, `e2^2`.
- Verified numerically in WSL: `tests/framework/test_otilib_spring_sensitivity.py`.

**visual:** three-panel strip — the 4-line residual, the 8-line `resasm.yml`, and the
resulting `public/summary.md` parameter-ranking table side by side.

---

## Slide 10 — Example: private black-box solver

- `templates/user_blackbox_residual/` — the model lives entirely inside the user's
  executable. `resasm.yml` says only:
  `residual: {type: executable, command: python my_solver.py --request {request} --response {response}}`
  and `tangent: {type: response}`.
- Per order `p`, the framework writes `request.json` (`u`, `parameters`,
  `direction_map`, `u_star_coefficients`, ...) and reads back `response.npz` with
  `R_order_<p>` and, optionally, `tangent`.
- The framework does the linear solve; the executable does the differentiation. For
  order ≥ 2 the executable must consume `u_star_coefficients` (the already-solved
  lower-order `U^(k)`) to rebuild `u*`.
- **OTILib is not needed on our side for this path.** We never link to, compile, or
  read the user's code.
- Honest reporting: the black-box never returns the real residual, so the run reports
  `residual_free_norm: null` with the reason, and the summary states that equilibrium
  was not verified.

**visual:** a sequence diagram — resasm → `request.json` → user's executable (drawn as
an opaque box) → `response.npz` → resasm → solve → package; the executable's interior
shaded solid black.

---

## Slide 11 — Current status

- **The user-facing workflow exists and runs**: `resasm init / check / run / report`
  (`residual_core/ui/cli.py`), four templates (python, blackbox, cpp, fortran),
  private/public output split.
- **OTILib**: external, GPLv3, **not vendored** — install with
  `scripts/setup_otilib.sh`. It does **not** build natively on Windows.
- **Verified**: on WSL, `bash scripts/run_otilib_tests_wsl.sh` runs the three OTILib
  test files with `RUN_OTILIB_TESTS=1` (so a skip would fail) — **6 passed**.
- **Known limitation**: on Windows without WSL the Python residual path cannot run at
  all — `oti_global.solve_python` hard-requires OTILib, so the advertised `dual1`
  backend is unreachable there. The black-box / C++ / Fortran paths are unaffected.
- **The FE framework (`residual_core/`)** is a separate, formulation-agnostic residual
  assembler. Its CP C3D8+UMAT backend is **verified offline** to machine precision
  (patch tests 1e-13, frame objectivity 5e-16, FD tangent 1.6e-16) — but the
  **Abaqus comparison is still pending**: Abaqus is not installed here, so ODB export,
  real-UMAT replay (needs ifort), and the reaction-force comparison are built and
  ready to run, **not executed**. They are marked pending, never passed.
- Crystal plasticity is one example backend, not the organizing principle.

**visual:** a status table — rows: user workflow, OTILib engine, Python path on
Windows, FE assembler core, CP C3D8+UMAT offline, CP vs Abaqus; columns: "verified /
pending / not supported", with the evidence (test file or script) in a third column.

---

## Slide 12 — What is next

- Close the Abaqus loop: run `Compression111` with the UMAT, export `fields.json`,
  confirm the C3D8 integration-point ordering, compare the assembled residual to `RF`.
- Build the real UMAT with ifort + Abaqus (MKL) and run the increment-by-increment
  replay against `S` / `SDV`.
- Implement the reserved `validation.solution_finite_difference_solver`: a true
  solution-level FD check that re-runs the user's nonlinear solver at perturbed
  parameters — the one check that can catch a wrong tangent.
- Make the Windows story better: today the Python path requires WSL.
- Timing study: one OTI evaluation per order versus one nonlinear re-solve per
  parameter for finite differences.
- Wire the OTI path into the FE element backends (today the FE assembler and the OTI
  global-residual path are separate).

**visual:** a two-column roadmap — "closes an open verification gap" (Abaqus, ifort,
solution-level FD) versus "extends reach" (Windows, timing study, FE + OTI
integration) — with the current release marked as a milestone on the left edge.
