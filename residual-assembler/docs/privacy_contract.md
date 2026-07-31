# Privacy contract

Every run writes one output directory (`output.dir`, default `resasm_output/`)
split in two:

```
resasm_output/
    private/     full numerical arrays — keep local
    public/      norms, rankings, status, timing — safe to share
```

`resasm_user` performs **no network I/O**: nothing is uploaded, and there is no
telemetry. "Private" therefore means *never written into `public/`* — the directory
you are expected to hand to a collaborator. Everything else stays on your disk
under your control.

Written by [`resasm_user/output.py`](../resasm_user/output.py).

## Private by default

None of these are written to `public/`:

| Thing | Where it lives |
|---|---|
| Your source code, mesh, geometry, material law | only in your own files — the framework never copies them anywhere |
| The full residual vector `R(u, a)` | `private/residual_real.npz` (key `residual`) — written **only when the residual was actually available**; on the black-box path it is omitted and replaced by `private/residual_real_UNAVAILABLE.json` rather than a fake zero vector |
| The tangent matrix `T` | `private/tangent.npz` (key `tangent`) — written only if a tangent was obtained |
| The OTI residual coefficients `R^(p)` and RHS `-R^(p)` | `private/rhs_order<p>.npz` (keys `residual` = `R^(p)`, `rhs` = `-R^(p)`) |
| The solution sensitivities `U^(p)` | `private/solution_sensitivities_order<p>.npz` (key `U`) |
| State / history variables | only in your `state.file` / `state.values` — never written to any output file |
| DOF partition (free / prescribed indices) | `private/dof_map.json` |
| Parameter → basis index map | `private/parameter_map.json` |
| Full validation record (norms, FD detail, backend diagnostics) | `private/validation_full.json` |
| Run metadata (name, order, backend, tangent source, parameter names) | `private/metadata.json` |

## Shareable by default — `resasm_output/public/`

Exactly five files, and no others:

| File | Contents |
|---|---|
| `summary.md` | Human-readable overview: run name, DOF count, parameter **names**, derivative order, backend, tangent source, residual free-DOF norm (or `n/a`), orders solved, the order-1 parameter ranking, and the RHS FD cross-check result if it ran. |
| `validation_summary.json` | Machine-readable: `residual_free_norm` (or `null`), `residual_free_norm_available`, `tangent_source`, `orders_solved`, `parameters` (names), `status`, `rhs_finite_difference_check`, `solution_finite_difference_check` (always `null` — not implemented), plus `residual_free_norm_reason` / `notes` when present. |
| `sensitivity_norms.csv` | One row per (order, direction): `order, direction, solution_sensitivity_norm, rhs_norm`. Norms only — never the vectors. |
| `parameter_ranking.csv` | `rank, parameter, order1_sensitivity_norm` — parameters sorted by ‖U^(1) column‖. |
| `timing.json` | `total_seconds`, `orders`, `unknowns`. |

## The nuance that matters: names yes, values no

**Parameter names are public. Parameter values are not.**

Verified against `resasm_user/output.py`:

- `public/summary.md` prints the parameter **names** (`| parameters | C11, C12, C44 |`)
  and their order-1 sensitivity **norms** in the ranking table.
- `public/parameter_ranking.csv` holds `parameter` (name) and
  `order1_sensitivity_norm` (a float).
- `public/validation_summary.json` holds `parameters` (the name list).
- `public/sensitivity_norms.csv` holds direction labels (`e1`, `e1*e2`, `e2^2`, …),
  which encode the *index* of a parameter — never its value.

The numeric parameter values you wrote in `resasm.yml` are **not written to any
output file — not even to `private/`**. `private/metadata.json` records only the
name list; `private/parameter_map.json` records only `name -> basis index`. Your
parameter values exist in exactly two places: your own `resasm.yml`, and the
transient `request.json` handed to your own black-box executable (see below).

## What `public/` does reveal (be honest with yourself)

`public/` is a lossy summary, not zero-knowledge. Before sharing it, know that it
discloses: the number of DOFs, the number and names of your parameters, the
derivative order and backend used, the tangent source, the residual norm on the
free DOFs (a scalar), and the magnitude (norm) and relative ranking of each
sensitivity direction. If a *parameter name* is itself confidential, rename it in
`resasm.yml` before running.

## Black-box mode — the strongest privacy path

With `residual.type: executable` the framework never imports, reads, or sees your
code. The entire exchange is two files in a temporary directory that is deleted
when the call returns
([`resasm_user/providers.py`](../resasm_user/providers.py)):

- **in** — `request.json`: `u`, `parameters`, `seed_directions`, `basis_count`,
  `truncation_order`, `order`, `direction_map`, `u_star_coefficients`, `state`,
  `time`, `dtime`. This is written *by* the framework *for* your solver, on your
  machine, and never leaves it.
- **out** — `response.npz` (or `.json`): `R_order_<p>` arrays, an optional
  `tangent`, and optional `diagnostics`.

The framework then only performs the linear solve `T U^(p) = -R^(p)`. Your model,
mesh, material law, and residual implementation stay entirely inside your
executable.

### The honest cost: no equilibrium verification

A black-box solver returns only the perturbed RHS coefficients. It **never returns
the real residual `R(u, a)`**, so there is nothing to take a norm of. The code
therefore records `residual_free_norm: null` with the reason
`"black-box did not expose real residual"`, and the public report says so instead
of printing a fabricated `0.000e+00`:

- `public/summary.md`: `| residual free-DOF norm | n/a — black-box did not expose real residual |`
- plus an explicit block: **"Equilibrium was NOT verified by this run."** … *"Nothing
  here should be read as confirming the supplied solution is converged."*
- `public/validation_summary.json`: `"residual_free_norm": null`,
  `"residual_free_norm_available": false`, `"residual_free_norm_reason": "..."`.

Also note: the RHS finite-difference cross-check does **not** run in black-box mode
(it needs to re-evaluate your residual at perturbed parameters). And no
solution-level finite-difference check exists on any path.

## How to audit `public/` before you share it

1. **Count the files.** There must be exactly five, and nothing else.
   ```bash
   ls resasm_output/public/
   # parameter_ranking.csv  sensitivity_norms.csv  summary.md  timing.json  validation_summary.json
   ```
2. **No arrays.** There must be no `.npz`, `.npy`, mesh, or source files.
   ```bash
   find resasm_output/public -type f ! -name '*.csv' ! -name '*.json' ! -name '*.md'
   ```
   (PowerShell: `Get-ChildItem resasm_output\public -Recurse -Exclude *.csv,*.json,*.md`)
3. **No parameter values.** Grep the whole directory for each numeric value you put
   under `parameters:` in `resasm.yml`. Nothing should match.
   ```bash
   grep -r "168400" resasm_output/public/    # expect: no output
   ```
4. **Read `summary.md` end to end.** It is short. Confirm the only numbers in it are
   DOF count, order, norms, and the ranking.
5. **Check the honesty flags** in `validation_summary.json`: if
   `residual_free_norm_available` is `false`, do not let a reader infer that
   equilibrium was verified — the `n/a` note in `summary.md` says so, keep it.
6. **Check `notes`.** If you set `validation.solution_finite_difference_solver`, a
   note stating that it is reserved and NOT implemented is added; keep it in the
   shared copy.
7. **Never share `private/`** unless you have deliberately decided to. `summary.md`
   only refers to it by relative path (`../private/`); it does not embed it.

See also: [simple_config_contract.md](simple_config_contract.md),
[glossary.md](glossary.md), [user_output_contract.md](user_output_contract.md).
