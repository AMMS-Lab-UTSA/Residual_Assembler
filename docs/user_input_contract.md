# User input contract

The **minimum** the framework needs to compute sensitivities.

## Required for every run

1. **Converged solution vector `u`** — a length-`unknowns` array (`solution.file`
   as `.npy`, or `solution.values`).
2. **Parameter names and real values** — `parameters: {name: value}`. These are
   the design variables that get seeded (`a_i* = a_i + e_i`).
3. **A residual evaluator** — one of:
   - Python: `residual(u, params, state=None, time=None) -> R`
   - Executable: a command that reads `request.json`, writes `response.npz`
   - Element: a Python module the framework assembles element-by-element
4. **A tangent** — `T = dR/du` at `u`. Provide it as a Python `tangent(...)`
   function, a saved `tangent.npz` file, or returned in the black-box response.
5. **Requested derivative order** — `sensitivity.order` (≥ 1).

## Optional

- `state` / `time` — history and time passed through to your residual.
- `constraints.free` / `constraints.prescribed` — DOF partition (default: all
  free; prescribed DOFs get zero sensitivity).
- `validation.rhs_finite_difference_check` — cross-check the **RHS / residual
  derivative** (`d(residual)/d(parameter)` at fixed `u`, solved with the same
  tangent). It is *not* a solution-level finite-difference validation and does not
  re-solve the nonlinear problem. (`validation.finite_difference` is a deprecated
  alias; `validation.solution_finite_difference_solver` is reserved and not
  implemented.)
- `output.dir` — where results go (default `resasm_output/`).
- parameter scaling — scale values in your residual if needed.

## You do NOT need mesh/element/material details

> The framework does not need to know your mesh if your residual function returns
> the **global** residual. It only needs mesh/element data if you want the
> framework to assemble the residual element-by-element (`residual.type: element`).

This is what lets private-code users keep everything local: expose only a global
`R(u, params)` (Python) or an executable that returns `R^(p)` (black-box).

## Black-box request/response

The executable is called with two file paths (`{request}`, `{response}`).

**request.json** (framework → executable):

```json
{ "schema": "resasm-user-request/1",
  "u": [ ... ],
  "parameters": { "k": 2.0, "f": 16.0 },
  "seed_directions": { "k": 1, "f": 2 },
  "basis_count": 2,
  "truncation_order": 1,
  "order": 1,
  "direction_map": { "1": [[1,0],[0,1]] },
  "u_star_coefficients": {},
  "state": null, "time": [0,0], "dtime": 0 }
```

**response.npz** (executable → framework): arrays `R_order_<p>` of shape
`(unknowns, N^(p))`, optional `tangent` `(unknowns, unknowns)`, optional
`diagnostics` (a JSON string). A `response.json` with
`residual_coefficients_by_order` + `tangent` + `diagnostics` is also accepted.

For `order ≥ 2`, also consume `u_star_coefficients` (the solved lower-order
`U^(k)`) to rebuild `u*` before differentiating.
