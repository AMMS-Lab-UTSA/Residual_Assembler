# Glossary

Terms as this repository actually uses them, with the concrete config key or
output artifact each one maps to.

## The problem

**residual** — The vector-valued function `R(u, a)` whose root is your converged
solution: `R = 0` at equilibrium. You supply it, and it is the *only* thing you
must expose. Configured by the `residual` block (`type`, `module` + `function`, or
`command`). Its real value at the converged `u` is stored in
`private/residual_real.npz` (key `residual`) — but only when the provider actually
returns it. A black-box provider never does, so on that path the file is omitted
and `private/residual_real_UNAVAILABLE.json` explains why (no zero vector is
written, because that would falsely read as verified equilibrium).

**converged solution** — The real solution vector `u` that already satisfies
`R(u, a) = 0` for your current parameter values. The toolkit never solves the
nonlinear problem; you hand it a converged `u` via `solution.file` or
`solution.values`, and its length defines `problem.unknowns` when that key is
omitted.

**tangent** — The Jacobian `T = dR/du` evaluated at the converged real solution.
It is the matrix on the left of `T U^(p) = -R^(p)` and is reused, unchanged, for
every order `p`. Configured by the `tangent` block (`python`, `file`, or
`response`) and stored in `private/tangent.npz` (key `tangent`).

**parameter** — A real design/material quantity you want sensitivities with
respect to (e.g. `C11`, `k`, `f`). Listed as `name: value` under `parameters`; the
key order fixes each parameter's imaginary basis index. Names are public; values
are never written to any output file.

**history-dependent model** — A residual that also depends on state/history
variables and on time, e.g. a plasticity law. Supplied through `state.file` /
`state.values` and `time.time` / `time.dtime`, which are passed to your residual
(and into the black-box `request.json`). The toolkit evaluates at that one fixed
state and time — it does not integrate a load path.

## The hypercomplex algebra

**OTI basis** — Order-Truncated-Imaginary algebra assigns one imaginary basis
element `e_i` to each parameter `a_i` (1-based, in `parameters` key order). Every
parameter is seeded simultaneously as `a_i* = a_i + e_i`, so a single OTI
evaluation of the residual carries all mixed derivatives at once. The mapping is
recorded in `private/parameter_map.json` and sent to a black-box solver as the
`seed_directions` field of the request.

**truncation order** — The maximum derivative order `q` the algebra retains;
everything above `q` is truncated away. Set by `sensitivity.order`; orders
`p = 1 … q` are generated and solved in sequence. Also called `nt` inside
`residual_core/algebra/otilib_adapter.py`.

**direction map** — The list of order-`p` directions, each an exponent multi-index
`κ` of length `m` (number of parameters). For `m = 2`, order 2 gives
`[2,0], [1,1], [0,2]`, labelled `e1^2`, `e1*e2`, `e2^2`. It is the `direction_map`
field of the black-box `request.json` (`{"p": [[exponents], ...]}`), and its labels
are the `direction` column of `public/sensitivity_norms.csv`. Each direction also
has a *recovery factor* `∏ᵢ κᵢ!` that converts a raw OTI coefficient into a true
partial derivative (it is 1 for every order-1 direction).

**R^(p)** — The order-`p` OTI coefficient matrix of the residual: shape
`(ndof, number of order-p directions)`, one column per direction. It is what the
tool extracts from the OTI residual (Python path) or asks the black-box to return
(`R_order_<p>` in the response). Stored in `private/rhs_order<p>.npz` under key
`residual`.

**rhs^(p)** — Simply `-R^(p)`: the right-hand side of the linear system
`T U^(p) = -R^(p)`. Stored in the same file, `private/rhs_order<p>.npz`, under key
`rhs`. Its per-direction norms are the `rhs_norm` column of
`public/sensitivity_norms.csv`.

**U^(p)** — The order-`p` solution sensitivity coefficient matrix, obtained by
solving `T U^(p) = -R^(p)` on the free DOFs (prescribed rows are zero). Stored in
`private/solution_sensitivities_order<p>.npz`, which ships **both** conventions:
`U_coefficients` (raw OTI Taylor coefficients; also the legacy key `U`) and
`U_derivatives` (= `recovery_factor × coefficient`, the true partial derivatives),
alongside `recovery_factors` and `direction_exponents`. At order 1 the factor is 1,
so `U^(1)[:, i]` *is* `du/da_i` — which is what `public/parameter_ranking.csv` ranks
by norm. Above order 1, use `U_derivatives`.

## The providers

**residual provider** — The adapter that evaluates your residual for the
framework. `residual.type: python` (and `element`, which currently uses the same
Python provider) loads `residual.module` and calls
`residual(u, params, state, time)`; the framework runs the OTI algebra itself. See
`resasm_user/providers.py::PythonResidual`.

**black-box provider** — `residual.type: executable`. The framework writes
`request.json`, runs your `residual.command`, and reads back `response.npz`
(`R_order_<p>` arrays, optional `tangent`, optional `diagnostics`). Your solver
generates `R^(p)` itself; the framework only does the linear solve. Your source is
never read. See `resasm_user/providers.py::ExecutableResidual`.

**tangent provider** — Where `T` comes from. `tangent.type: python` calls a
`tangent(u, params, state, time)` function in your residual module;
`tangent.type: file` loads an `(n, n)` array from `.npz`/`.npy`;
`tangent.type: response` takes the `tangent` array out of the black-box response.
The choice is echoed as `tangent_source` in `public/validation_summary.json` and
`public/summary.md`.

## Validation

**RHS finite-difference check** — **Implemented.** Enabled with
`validation.rhs_finite_difference_check: true` (legacy alias:
`validation.finite_difference`), Python/`element` path only. It central-differences
`d(residual)/d(parameter)` **at the fixed solution `u`**, solves with the **same**
tangent the hypercomplex path used, and compares to `U^(1)`. It therefore validates
the generated RHS and that solve — it cannot detect an error in the tangent itself
(a wrong `T` sits on both sides and cancels). Results appear as
`rhs_finite_difference_check` in `public/validation_summary.json` and as a table in
`summary.md`.

**solution finite-difference check** — **NOT implemented.** This would be a
genuinely different thing: re-running *your nonlinear solver* at perturbed
parameters `a_i ± h` and finite-differencing the resulting converged solutions —
the only check that can catch a wrong tangent. `validation.solution_finite_difference_solver`
is a **reserved config field only**; setting it performs no check and instead adds
an explicit note saying so. `public/validation_summary.json` always reports
`"solution_finite_difference_check": null`. Do not confuse it with the RHS check
above.

## The outputs

**private package** — `resasm_output/private/`. Full numerical arrays, kept local:
`residual_real.npz`, `tangent.npz`, `rhs_order<p>.npz`,
`solution_sensitivities_order<p>.npz`, plus `metadata.json`, `parameter_map.json`,
`dof_map.json`, and `validation_full.json`. This is the answer you actually use.

**public report** — `resasm_output/public/`. Exactly five files, safe to share:
`summary.md`, `validation_summary.json`, `sensitivity_norms.csv`,
`parameter_ranking.csv`, `timing.json`. Norms, rankings, status, and timing only —
no source, mesh, full residual, tangent, state, or parameter values. See
[privacy_contract.md](privacy_contract.md).

**recovery factor** — `Π_i (κ_i!)` for a direction with exponents `κ`. Converts a raw
OTI Taylor **coefficient** into a true partial **derivative**:
`derivative = recovery_factor × coefficient`. It is 1 at order 1 and for mixed
directions, and >1 for repeated ones (`d2/dk2` → 2, `d3/dk3` → 6). Exported per
column in `private/direction_map_order<p>.json` and as `recovery_factors` inside
each order's `.npz`.

**OTI coefficient vs recovered derivative** — Two different numbers above order 1.
`U_coefficients` / `residual_coefficients` are Taylor coefficients (what OTI emits);
`U_derivatives` / `residual_derivatives` are the partial derivatives you almost always
want. Public reports quote **derivatives**; raw coefficients are kept for traceability.
The legacy keys `U`, `residual`, `rhs` hold **coefficients**.
