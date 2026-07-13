# Input objects

Everything the toolkit needs to build the **residual-sensitivity execution
contract** — the ingredients of

```
T U^(p) = -R^(p)
```

where `T` is the tangent at your converged real solution, `R^(p)` are the
order-`p` OTI coefficients of your residual, and `U^(p)` are the order-`p`
solution sensitivities.

All inputs come from one file, `resasm.yml` (schema:
[`schemas/resasm_config_v1.schema.json`](../schemas/resasm_config_v1.schema.json),
loader: `resasm_user/config.py`). Paths inside it are resolved relative to the
directory containing the config.

---

## Two things to know before you start

> **`problem.unknowns` is OPTIONAL.** It is inferred from the solution vector
> (`solution.file` / `solution.values`). If you give both, they must agree —
> a mismatch is a hard error, never a silent truncation
> (`config.py::_resolve_unknowns`).

> **You do NOT need to provide mesh, element, or material details** if your
> residual provider already returns the **global** residual `R(u, params)`.
> The framework never asks for a mesh on the `python` or `executable` paths. It
> only needs the vector `u`, the parameter values, and something that can
> evaluate `R` (and something that can supply `T`).

---

## Object table

| Object | Where it comes from | Shape / type | Required? | If missing (actual error) |
|---|---|---|---|---|
| Job name | `problem.name` | string | yes | `Missing required field: problem.name` |
| Number of unknowns | `problem.unknowns` | int ≥ 1 | **no** — inferred | Inferred from the solution vector. If neither is determinable: `cannot determine the number of unknowns: the solution vector is not readable yet and problem.unknowns is not set` |
| Solution vector `u` | `solution.file` (`.npy`) or `solution.values` | `(ndof,)` float, flattened with `.ravel()` | yes | `Missing required field: solution` / `solution needs a converged solution vector` / at run time `solution file not found: <path>` |
| Parameter map | `parameters: {name: real}` | mapping, ≥ 1 entry, values coerced with `float()` | yes | `Missing required field: parameters` / `parameters must be a non-empty mapping of name -> real value` / `parameter 'k' must be a real number, got ...` |
| Residual provider | `residual.type` + `residual.module`/`function` or `residual.command` | see below | yes | `Missing required field: residual.type`; `residual.type 'x' is not supported`; `residual.module is required for a Python residual`; `residual.command is required for a black-box executable` |
| Tangent provider | `tangent.type` = `python` \| `file` \| `response` | `(ndof, ndof)` dense float | yes (some source must resolve) | Python path: `no tangent available for the Python path`. Black-box: `no tangent available: provide tangent.file or return \`tangent\` in the black-box response`. Wrong shape: `tangent has shape (a, b) but must be (ndof, ndof)` |
| Sensitivity request | `sensitivity.order`, `sensitivity.backend` | int ≥ 1; string | `order` yes, `backend` no (default `otilib`) | `Missing required field: sensitivity.order` / `sensitivity.order must be >= 1 (got 0)` |
| State / history | `state.file` (`.npy`, `allow_pickle=True`) or `state.values` | anything your residual understands | no | `state = None` is passed to your residual |
| Time / dtime | `time.time`, `time.dtime` | float, float | no | both default to `0.0` |
| DOF partition | `constraints.free` or `constraints.prescribed` | list of 0-based int indices | no | all DOFs free |
| Output directory | `output.dir` | string | no | defaults to `resasm_output` |
| Validation options | `validation.*` | see below | no | no cross-check is run |

---

## 1. Solution vector `u`

The converged **real** solution of your problem. This is the point at which the
tangent is taken and around which the sensitivities are expanded.

```yaml
solution:
  file: solution.npy        # a .npy array
# or:
solution:
  values: [2.0]             # an inline list
```

- Loaded by `runner.py::_load_solution`, then `np.asarray(..., float).ravel()`
  → shape `(ndof,)`. A 2-D `.npy` is flattened, not rejected.
- `ndof` (a.k.a. `problem.unknowns`) is defined by this vector.
- Nothing checks that `u` is actually converged. On the Python path the residual
  norm on the free DOFs is computed and reported (and `resasm check` warns if it
  exceeds `1e-4`). On the black-box path it cannot be checked at all — see
  [`output_objects.md`](output_objects.md).

## 2. Parameter map

```yaml
parameters:
  k: 2.0
  f: 16.0
```

These are the design variables that get seeded (`a_i* = a_i + e_i`).

**The order of the keys is the OTI basis / seed order.** Key #1 → imaginary
basis `e1`, key #2 → `e2`, and so on (`oti_global.py::solve_python`,
`providers.py::ExecutableResidual.eval_rhs` → `seed_directions`,
`output.py` → `parameter_map.json`). Every downstream index — the columns of
`R^(p)` / `U^(p)`, the exponent vectors in `direction_map`, the row order in
`parameter_ranking.csv`, and the `params[0]`/`params(1)` positions in the C++
and Fortran templates — follows this same order. Reordering the keys in
`resasm.yml` reorders all of it.

Values must be real scalars. Complex values, vectors, or per-element fields are
not supported as parameters.

## 3. Residual provider

`residual.type` accepts exactly three values (`config.py::_validate_residual`):

| `residual.type` | What it means | Required fields |
|---|---|---|
| `python` | A Python callable `residual(u, params, state=None, time=None) -> R`. The framework seeds OTI numbers and extracts `R^(p)` itself. | `residual.module` (path to a `.py`), `residual.function` (default `residual`) |
| `executable` | A local command that reads `request.json` and writes a response file containing `R^(p)`. Your model never leaves your process. C++ and Fortran providers use this path. | `residual.command`, containing the `{request}` and `{response}` placeholders |
| `element` | Accepted by the config, but in the current code it is handled **identically to `python`**: `runner.py` builds a `PythonResidual`, loading `residual.module` and calling `residual.function` (default `residual`). There is no mesh/element assembly in the run path. | `residual.module` |

```yaml
# Python
residual:
  type: python
  module: user_residual.py
  function: residual

# Executable (also used for compiled C++ / Fortran providers)
residual:
  type: executable
  command: ./my_solver --request {request} --response {response}
```

The full call/response semantics are in
[`residual_provider_contract.md`](residual_provider_contract.md).

## 4. Tangent provider

`T = dR/du` at the converged `u`, dense, shape `(ndof, ndof)`. This is what the
solve `T U^(p) = -R^(p)` uses (restricted to the free DOFs).

`runner.py::_acquire_tangent` acts on exactly two `tangent.type` values; the
others fall through to "no tangent loaded here":

| `tangent.type` | Behaviour |
|---|---|
| `python` | Calls `tangent(u, params, state, time)` in the **same module** as the residual (function name from `tangent.function`, default `tangent`). Only honoured when the provider is a `PythonResidual`. |
| `file` | `tangent.file` — a `.npz` (key `tangent`, else the first array in the archive) or a `.npy` (`providers.py::load_tangent_file`). |
| `response` | Nothing is loaded up front; the tangent must come back in the black-box response as the `tangent` array/field. If the `tangent` block is omitted entirely and `residual.type: executable`, the config defaults to `{"type": "response"}`. |
| `none` | Allowed by the schema, but inert: no tangent is loaded. On the Python path this ends as `no tangent available for the Python path`. |

```yaml
tangent:
  type: python
  function: tangent
# or
tangent:
  type: file
  file: tangent.npz
# or
tangent:
  type: response      # black-box returns it
```

The shape is checked (`_check_tangent`): anything other than `(ndof, ndof)` is a
`ConfigError`. The tangent is used as-is; it is never verified against the
residual (the RHS finite-difference check reuses the *same* `T` on both sides,
so it cannot detect a wrong tangent).

## 5. Sensitivity request

```yaml
sensitivity:
  order: 2
  backend: otilib
```

- `order` — the truncation order `q`. Orders `p = 1..q` are generated and solved.
- `backend` — recorded in `private/metadata.json` and printed by `resasm check`.
  **The run path always uses OTILib** for `python`/`element`
  (`oti_global.solve_python` imports `OtiContext` unconditionally and raises
  `OtiUnavailable` if OTILib is not importable). The schema also lists `dual1`,
  but no `dual1` code path exists in `resasm_user`: `resasm check` rejects it for
  `order > 1`, and `resasm run` ignores the field. Use `otilib`.
- The black-box path does not need OTILib on the framework side at all — the
  executable produces `R^(p)` and the framework only does the linear solve.

The number of directions per order is `N^(p) = C(p + m - 1, p)` for `m`
parameters, so `R^(p)` and `U^(p)` are `(ndof, N^(p))` matrices.

---

## Optional objects

### State / history

```yaml
state:
  file: state.npy         # np.load(..., allow_pickle=True)
# or
state:
  values: {...}           # any inline structure
```

Passed through to your residual untouched as the `state` argument (`None` if
absent). The framework never interprets it.

**Black-box caveat:** on the `executable` path the state is placed into
`request.json` via `json.dump`, so it must be JSON-serializable. A NumPy array
loaded from `state.file` is not, and the run will fail when serializing the
request. Use `state.values`, or keep the history inside your executable.

### Time

```yaml
time:
  time: 1.0
  dtime: 0.1
```

`runner.py::_time_pair` builds `((t, t), dtime)` from the single `time.time`
value. What your residual actually receives differs by call site — see
[`residual_provider_contract.md`](residual_provider_contract.md#the-time-argument).
`dtime` is only forwarded to the black-box request; the Python residual never
sees it.

### Free / prescribed DOF partition

```yaml
constraints:
  prescribed: [0, 7, 12]   # 0-based DOF indices
# or
constraints:
  free: [1, 2, 3, ...]
```

`runner.py::_free_mask`: if `prescribed` is present it wins (`free` is then
ignored); otherwise `free` lists the free DOFs. Default: **all DOFs free**.
The solve is done on the free–free block `T[free, free]`; prescribed DOFs get
zero sensitivity rows in every `U^(p)`.

### Validation

```yaml
validation:
  rhs_finite_difference_check: true
```

| Field | Status |
|---|---|
| `validation.rhs_finite_difference_check` | Implemented, **Python/element path only** (`runner.py::rhs_finite_difference_check`). Central-differences `d(residual)/d(parameter)` at the **fixed** `u`, solves with the **same** tangent, compares to `U^(1)`. It validates the generated RHS and that solve. It does **not** re-solve the nonlinear problem and cannot detect an error in `T`. |
| `validation.finite_difference` | Deprecated alias of the above, still accepted (`runner.py::_wants_rhs_fd`). Output always reports the accurate name. |
| `validation.solution_finite_difference_solver` | **RESERVED — NOT IMPLEMENTED.** A real solution-level FD check would re-run your nonlinear solver at perturbed parameters. Setting it performs no such check; the run adds an explicit note saying so (`runner.py::_solution_fd_note`) and `public/validation_summary.json` keeps `solution_finite_difference_check: null`. |
| `validation.fd_step` | Declared in the schema but **not read** by the current code. The FD step is hard-coded: `h = 1e-6 * max(1, |a_i|)`. |

### Output directory

```yaml
output:
  dir: resasm_output      # default
```

---

## Not supported (do not put these in `resasm.yml`)

- **Parameter scaling / normalisation.** There is no scaling field anywhere in
  `config.py` or `runner.py`. Parameter values are passed to your residual (and
  seeded) exactly as written. If you need scaling, do it inside your residual.
- **Sparse tangents.** `T` is loaded and inverted densely with
  `np.linalg.solve`.
- **Mesh / element / material blocks.** Not read on any implemented path.
- **Multiple load steps / increments.** One `(time, dtime)` pair, one run.
- **A `dual1` sensitivity backend.** Declared in the schema, absent from the code.

---

## Minimal complete example

```yaml
problem:
  name: spring_demo
  # unknowns omitted -> inferred from solution.npy

residual:
  type: python
  module: user_residual.py
  function: residual

tangent:
  type: python
  function: tangent

parameters:
  k: 2.0
  f: 16.0

solution:
  file: solution.npy

sensitivity:
  order: 2
  backend: otilib

validation:
  rhs_finite_difference_check: true
```

Then:

```
resasm check resasm.yml     # readiness report, stops at the first blocker
resasm run   resasm.yml     # writes resasm_output/private + resasm_output/public
```

Working, runnable versions of each path live in `templates/`
(`resasm init --template python|blackbox|cpp|fortran --out my_job`).
