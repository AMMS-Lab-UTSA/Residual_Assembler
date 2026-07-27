# Output objects

Every run writes one directory (default `resasm_output/`, override with
`output.dir`) split in two. Produced entirely by
`resasm_user/output.py::write_outputs`.

```
resasm_output/
    private/          full numerical arrays — keep local
        metadata.json
        parameter_map.json
        dof_map.json
        residual_real.npz
        tangent.npz
        rhs_order<p>.npz                    one per solved order
        solution_sensitivities_order<p>.npz one per solved order
        direction_map_order<p>.json         what each column means (+ recovery factors)
        validation_full.json
    public/           safe to share
        summary.md
        timing.json
        parameter_ranking.csv
        sensitivity_norms.csv
        validation_summary.json
```

---

## The hard rule

**`public/` must NEVER contain:**

- source code,
- mesh, element, or material data,
- the full residual vector,
- the full tangent matrix,
- internal state variables,
- private parameter **values**.

Parameter **names** and **norms** are the only model-derived quantities that
cross into `public/`. Everything with a shape bigger than a scalar stays in
`private/`. If you add anything to `public/`, it must satisfy this rule.

---

## Shapes and conventions used below

- `ndof` — number of unknowns.
- `m` — number of parameters.
- `N^(p) = C(p + m - 1, p)` — number of order-`p` directions, i.e. the number of
  columns of `R^(p)` and `U^(p)`.
- **Column order** for order `p` comes from
  `residual_core/algebra/otilib_adapter.py::order_directions(m, p)`:
  `combinations_with_replacement` over the basis indices `1..m`, where basis `i`
  is the `i`-th key of `parameters`. Labels look like `e1`, `e2`, `e1*e2`,
  `e1^2`.
- **Coefficient convention.** `R^(p)` and `U^(p)` hold the **raw OTI
  coefficients** (`OtiContext.coeff` → the imaginary coefficient), not the
  partial derivatives. The true partial derivative along a direction is
  `coefficient × recovery_factor`, where `recovery_factor` is the product of the
  factorials of that direction's exponents (`e1^2` → `2! = 2`, `e1*e2` → `1`).
  At order 1 every factor is 1, so `U^(1)[:, i] = du/da_i` exactly.
- Prescribed DOFs (`constraints.prescribed`) are **zero rows** in every `U^(p)`;
  the solve only touches the free–free block of `T`.

---

## private/

Full arrays. This is the part that reveals your model — it never leaves your
machine unless you send it.

### `metadata.json`

Run identity, so any later reader knows what produced the arrays.

| Key | Value |
|---|---|
| `name` | `problem.name` |
| `unknowns` | `ndof` (int) |
| `order` | requested `sensitivity.order` |
| `backend` | `sensitivity.backend` as written in the config |
| `residual_type` | `python` \| `executable` \| `element` |
| `tangent_source` | `python` \| `file` \| `response` \| `none` — where `T` actually came from |
| `parameters` | list of parameter names, **in seed order** (names only, no values) |
| `created` | local timestamp, `%Y-%m-%dT%H:%M:%S` |
| `engine` | `"otilib"` (Python path) or `"executable"` (black-box), or `null` |

### `parameter_map.json`

Parameter name → 1-based OTI imaginary basis index.

```json
{"k": 1, "f": 2}
```

This is the key to reading every column index and every exponent vector.

### `dof_map.json`

```json
{"ndof": 1, "free_dofs": [0], "prescribed_dofs": []}
```

0-based DOF indices, derived from `constraints`.

### `residual_real.npz` — *written only when the residual actually exists*

| Key | Shape | Content |
|---|---|---|
| `residual` | `(ndof,)` | `R(u, params, state, time)` at the converged real solution |

**Black-box:** on the `executable` path the provider never hands back the real
residual, so **this file is not written at all**. In its place you get:

### `residual_real_UNAVAILABLE.json` (black-box only)

| Key | Content |
|---|---|
| `residual_available` | `false` |
| `reason` | e.g. `"black-box did not expose real residual"` |
| `note` | why no zero vector is written |

The tool deliberately does **not** write a zero vector here. A `residual_real.npz`
full of zeros would read as `‖R‖ = 0` — i.e. as if equilibrium had been verified —
when in fact the residual was never returned. The private package refuses that lie
for the same reason the public report prints `n/a` instead of `0.000e+00`.

### `tangent.npz`

| Key | Shape | Content |
|---|---|---|
| `tangent` | `(ndof, ndof)` | the tangent `T = dR/du` actually used in the solve — whether it came from your Python `tangent()`, a file, or the black-box response |

Written only when a tangent is available (a run without one raises before this
point, so in practice it is always present).

### ⚠ Coefficients vs derivatives — read this before using any order ≥ 2 array

The OTI evaluation produces **Taylor coefficients**, not partial derivatives. They
are related by the **recovery factor**:

```
    derivative = recovery_factor * coefficient ,    recovery_factor = Π_i (κ_i !)
```

It is `1` for every order-1 direction and for mixed directions (all `κ_i = 1`), and
greater than 1 whenever a parameter repeats — e.g. `d2/dk2` → `2! = 2`,
`d3/dk3` → `3! = 6`.

**You never have to apply it yourself.** Every order exports *both* conventions, plus
the map that connects them. `*_coefficients` are kept for traceability;
`*_derivatives` are what you almost certainly want; **the public report quotes
derivatives.**

### `rhs_order<p>.npz` — one per solved order `p = 1..order`

| Key | Shape | Content |
|---|---|---|
| `residual_coefficients` | `(ndof, N^(p))` | `R^(p)` — raw order-`p` Taylor coefficients |
| `rhs_coefficients` | `(ndof, N^(p))` | `-R^(p)` — the right-hand side of `T U^(p) = -R^(p)` |
| `residual_derivatives` | `(ndof, N^(p))` | `recovery_factor · R^(p)` |
| `rhs_derivatives` | `(ndof, N^(p))` | `-recovery_factor · R^(p)` |
| `recovery_factors` | `(N^(p),)` | `Π κ_i!` per column |
| `direction_exponents` | `(N^(p), m)` | the exponent vector `κ` per column |
| `residual` *(legacy)* | `(ndof, N^(p))` | **raw coefficients** — same as `residual_coefficients` |
| `rhs` *(legacy)* | `(ndof, N^(p))` | **raw coefficients** — same as `rhs_coefficients` |

> The legacy keys `residual` / `rhs` are preserved for backward compatibility and
> hold **raw coefficients**, not derivatives. New code should read the explicit
> `*_coefficients` / `*_derivatives` keys.

### `solution_sensitivities_order<p>.npz` — one per solved order

| Key | Shape | Content |
|---|---|---|
| `U_coefficients` | `(ndof, N^(p))` | `U^(p)` — raw Taylor coefficients; solution of `T_ff U^(p)_f = -R^(p)_f`, zero on prescribed DOFs |
| `U_derivatives` | `(ndof, N^(p))` | `recovery_factor · U^(p)` — **the actual partial derivatives** |
| `recovery_factors` | `(N^(p),)` | `Π κ_i!` per column |
| `direction_exponents` | `(N^(p), m)` | the exponent vector `κ` per column |
| `U` *(legacy)* | `(ndof, N^(p))` | **raw coefficients** — same as `U_coefficients` |

At **order 1** every factor is 1, so `U == U_coefficients == U_derivatives`, and
column `i` is `du/da_i`. At **order 2**, e.g. for the spring `R = k u³ − f` at
`k=2, f=16, u=2`:

| column | factor | `U_coefficients` | `U_derivatives` | exact |
|---|---|---|---|---|
| `d2/dk2` | 2 | `+0.111111111` (= 1/9) | `+0.222222222` | `d²u/dk² = 2/9` |
| `d2/dk_df` | 1 | `−0.006944444` | `−0.006944444` | (factor 1) |
| `d2/df2` | 2 | `−0.000868056` | `−0.001736111` | `d²u/df² = −1/576` |

Reading `U` (or `U_coefficients`) as a second derivative would under-report `d2/dk2`
by a factor of 2. Pinned by `tests/framework/test_oti_recovery_factor.py`.

### `direction_map_order<p>.json` — one per solved order

Tells you what each **column** means. No column ordering has to be inferred.

| Key | Content |
|---|---|
| `order` | `p` |
| `parameter_names` | the parameter list, in basis order |
| `note` | the `derivative = recovery_factor * coefficient` rule |
| `columns[]` | `index`, `exponents`, `label` (e.g. `d2/dk2`), `oti_label` (e.g. `e1^2`), `recovery_factor`, `parameter_names` |

```json
{
  "order": 2,
  "parameter_names": ["k", "f"],
  "columns": [
    {"index": 0, "exponents": [2, 0], "label": "d2/dk2",    "recovery_factor": 2.0},
    {"index": 1, "exponents": [1, 1], "label": "d2/dk_df",  "recovery_factor": 1.0},
    {"index": 2, "exponents": [0, 2], "label": "d2/df2",    "recovery_factor": 2.0}
  ]
}
```

### `validation_full.json`

The complete `validation` dict assembled by `runner.py`:

| Key | When | Content |
|---|---|---|
| `status` | always | `"ok"` |
| `residual_free_norm` | always | `‖R(u)[free]‖` on the Python/element path; `null` on the black-box path |
| `residual_free_norm_reason` | black-box only | `"black-box did not expose real residual"` |
| `rhs_finite_difference_check` | when `validation.rhs_finite_difference_check` is on (Python path) | `{max_rel_error, status: "pass"|"warn", checks: "..."}` — `pass` iff `max_rel_error < 1e-4` |
| `notes` | when `validation.solution_finite_difference_solver` was set | a list containing the explicit "reserved but NOT implemented" note |

Note: the diagnostics your black-box returns are parsed
(`providers.py::_read_response`) and collected in memory, but only the
`engine` string reaches disk (in `metadata.json`). The rest of the executable's
`diagnostics` payload is not persisted to any output file.

---

## public/

Norms, rankings, status, timing. For a collaborator, a reviewer, or a paper —
nothing here reconstructs your model.

### `summary.md`

Human-readable report. Contains:

1. **Run table** — unknowns, parameter *names*, derivative order, backend,
   tangent source, residual free-DOF norm, orders solved.
2. **Parameter ranking table** — rank, parameter name, `‖du/da‖`.
3. **RHS finite-difference cross-check** section, *only if* the check ran: max
   relative error and pass/warn, with an explicit paragraph saying it checks
   `d(residual)/d(parameter)` at the fixed `u` with the same tangent, does not
   re-solve the nonlinear problem, and therefore cannot detect an error in the
   tangent itself.
4. Any `notes` (e.g. the reserved-and-not-implemented solution FD note).
5. A file list.

**Black-box honesty.** When the residual was never exposed, the residual-norm
cell is *not* `0.000e+00`. It prints

```
| residual free-DOF norm | n/a — black-box did not expose real residual |
```

followed by:

> **Equilibrium was NOT verified by this run.** The residual vector was never
> available (black-box did not expose real residual), so no residual norm could
> be computed. Nothing here should be read as confirming the supplied solution
> is converged.

### `timing.json`

```json
{"total_seconds": 0.0421, "orders": 2, "unknowns": 1}
```

Wall-clock for the whole `run_from_config` call, the requested order, and `ndof`.

### `parameter_ranking.csv`

Which parameters matter, by order-1 solution-sensitivity magnitude.

| Column | Content |
|---|---|
| `rank` | 1 = most influential |
| `parameter` | parameter **name** (no value) |
| `order1_sensitivity_norm` | `‖U^(1)[:, i]‖₂`, `%.6e` |

Sorted by norm, descending. If order 1 was not solved, the norms are `0.0`.

### `sensitivity_norms.csv`

One row per (order, direction).

| Column | Content |
|---|---|
| `order` | `p` |
| `direction` | human-readable derivative label — `d/dk`, `d2/dk2`, `d2/dk_df`, … |
| `oti_direction` | the raw OTI label — `e1`, `e1^2`, `e1*e2`, … |
| `recovery_factor` | `Π κ_i!` used for this column (1 at order 1) |
| `solution_sensitivity_norm` | `‖U^(p)_derivatives[:, col]‖₂`, `%.6e` — **recovered derivative**, not raw coefficient |
| `rhs_norm` | `‖R^(p)_derivatives[:, col]‖₂`, `%.6e` — likewise recovered |

> **Public reports quote recovered derivatives.** Order-1 rows are numerically
> identical to the old raw-coefficient convention; order ≥ 2 rows for repeated
> directions are larger by their factorial, which is the correct value.

Norms only. The vectors themselves stay in `private/`.

### `validation_summary.json`

Machine-readable status (this is what `resasm report` and
`resasm_user.read_report` read).

| Key | Content |
|---|---|
| `residual_free_norm` | float, **or `null`** when the residual was never available |
| `residual_free_norm_available` | `true` / `false` — `false` means no equilibrium check was possible |
| `residual_free_norm_reason` | present only when the norm is `null`, e.g. `"black-box did not expose real residual"` |
| `tangent_source` | `python` \| `file` \| `response` \| `none` |
| `orders_solved` | e.g. `[1, 2]` |
| `parameters` | parameter names, in seed order |
| `status` | `"ok"` |
| `rhs_finite_difference_check` | the FD-check object, or `null` if it did not run |
| `solution_finite_difference_check` | **always `null`** — a solution-level FD check is not implemented |
| `notes` | present only when there is something to say (e.g. the reserved solution-FD field was set) |

**Black-box mode** produces exactly:

```json
{
  "residual_free_norm": null,
  "residual_free_norm_available": false,
  "residual_free_norm_reason": "black-box did not expose real residual",
  "tangent_source": "response",
  "orders_solved": [1],
  "parameters": ["k", "f"],
  "status": "ok",
  "rhs_finite_difference_check": null,
  "solution_finite_difference_check": null
}
```

A fabricated `0.000e+00` would imply the run verified equilibrium. It did not,
and it says so.

---

## Reading a run back

```python
from resasm_user import read_report
rep = read_report("resasm_output")
rep["validation_summary"]   # public/validation_summary.json
rep["metadata"]             # private/metadata.json
```

or `resasm report resasm_output/`.

To get the arrays:

```python
import numpy as np
z   = np.load("resasm_output/private/solution_sensitivities_order1.npz")
U1  = z["U_derivatives"]                                                           # (ndof, m) TRUE derivatives
#   z["U_coefficients"] are the raw Taylor coefficients; at order 1 they are equal.
#   For order >= 2 ALWAYS use U_derivatives (or multiply by z["recovery_factors"]).
R1  = np.load("resasm_output/private/rhs_order1.npz")["residual"]                  # R^(1)
b1  = np.load("resasm_output/private/rhs_order1.npz")["rhs"]                       # -R^(1)
T   = np.load("resasm_output/private/tangent.npz")["tangent"]                      # (ndof, ndof)
```
