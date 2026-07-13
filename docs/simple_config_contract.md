# Simple config contract (`resasm.yml`)

One small file describes the whole job. It is validated by
[`schemas/resasm_config_v1.schema.json`](../schemas/resasm_config_v1.schema.json)
and loaded by [`resasm_user/config.py`](../resasm_user/config.py).

Five blocks are required: `problem`, `residual`, `parameters`, `solution`,
`sensitivity`. Everything else is optional.

## Minimal Python residual

```yaml
problem:
  name: spring_demo
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
```

## Minimal black-box residual

```yaml
problem:
  name: private_solver_demo
residual:
  type: executable
  command: ./my_solver --request {request} --response {response}
tangent:
  type: response
parameters:
  C11: 168400.0
  C12: 121400.0
  C44: 75400.0
solution:
  file: converged_u.npy
sensitivity:
  order: 1
  backend: otilib
```

## `problem.unknowns` is optional

Neither example above declares `problem.unknowns`. It is **inferred from the
solution vector** (`solution.file` or `solution.values`). If you declare it *and*
a solution is readable, the two must agree — a mismatch is a **hard error**, never
a silent truncation:

```
problem.unknowns is 3 but the solution vector has 5 entries
```

If the solution vector cannot be read at load time *and* `problem.unknowns` is
absent, the load fails and asks you for one of the two.

## Paths

All relative paths (`residual.module`, `residual.command`, `solution.file`,
`tangent.file`, `state.file`, `output.dir`) are resolved against the directory
containing `resasm.yml`. The black-box command also runs with that directory as
its working directory.

## Key reference

| Key | Required? | Values | Meaning |
|---|---|---|---|
| `problem.name` | yes | string | Label for the run. Appears in `private/metadata.json` and `public/summary.md`. |
| `problem.unknowns` | no | integer ≥ 1 | Number of DOFs. Inferred from the solution vector when omitted; must agree with it when given. |
| `residual.type` | yes | `python` \| `executable` \| `element` | Which provider evaluates your residual. `element` is accepted but is currently handled by the **same** Python provider as `python` (see notes). |
| `residual.module` | for `python` / `element` | path to a `.py` | File containing your residual function. |
| `residual.function` | no | string (default `residual`) | Name of the residual callable in `residual.module`. Signature: `residual(u, params, state, time)` returning a length-`unknowns` vector. |
| `residual.command` | for `executable` | shell command | Must contain the `{request}` and `{response}` placeholders, e.g. `./my_solver --request {request} --response {response}`. |
| `tangent.type` | no | `python` \| `file` \| `response` \| `none` | Where the tangent `T = dR/du` comes from. Defaults to `{type: response}` when `residual.type: executable` and no `tangent` block is given. |
| `tangent.function` | for `tangent.type: python` | string (default `tangent`) | Tangent callable, looked up in `residual.module`. Signature `tangent(u, params, state, time)` returning an `(n, n)` array. |
| `tangent.file` | for `tangent.type: file` | path to `.npz` / `.npy` | For `.npz` the key `tangent` is used if present, otherwise the first array. Shape must be exactly `(n, n)`. |
| `parameters` | yes | non-empty map `name: real value` | The design parameters that get seeded. Order of keys fixes the imaginary basis index (first key → `e1`, …). |
| `solution.file` | one of the two | path to `.npy` | The converged real solution `u`. |
| `solution.values` | one of the two | list of numbers | Inline alternative to `solution.file` (takes precedence if both are present). |
| `sensitivity.order` | yes | integer ≥ 1 | Truncation order `q`: orders 1…`q` are generated and solved. |
| `sensitivity.backend` | no | `otilib` (default) \| `dual1` | Recorded in `private/metadata.json` and `public/summary.md`. See notes: the Python path always runs OTILib. |
| `constraints.prescribed` | no | list of integer DOF indices | Those DOFs are fixed; the rest are free. Takes precedence over `constraints.free`. |
| `constraints.free` | no | list of integer DOF indices | Only those DOFs are free; the rest are prescribed. Default (neither key): all DOFs free. Prescribed rows of every `U^(p)` are zero. |
| `state.file` | no | path to a `.npy` / `.npz` | History/state variables, loaded and passed to your residual as `state`. |
| `state.values` | no | any YAML value | Inline state, passed through unchanged. Used only if `state.file` is absent. |
| `time.time` | no | number (default `0.0`) | Passed through to your residual / the black-box request. |
| `time.dtime` | no | number (default `0.0`) | Time increment. Sent in the black-box request; **not** passed to a Python residual. |
| `output.dir` | no | path (default `resasm_output`) | Output root. Gets a `private/` and a `public/` subdirectory. |
| `validation.rhs_finite_difference_check` | no | bool | Run the RHS finite-difference cross-check (see below). Python/`element` path only. |
| `validation.finite_difference` | no | bool | **DEPRECATED** alias of `rhs_finite_difference_check`, same behaviour. The output always reports the accurate name. |
| `validation.solution_finite_difference_solver` | no | any | **RESERVED, NOT IMPLEMENTED.** Setting it performs no check; the run adds an explicit note saying so to `validation_summary.json` and `summary.md`. |
| `validation.fd_step` | no | number | Present in the schema but **not read** by `resasm_user`. The RHS check uses a fixed relative step of `1e-6`. |

## `validation.rhs_finite_difference_check` — what it actually checks

It central-differences the **residual** with respect to each parameter **at the
fixed solution `u`**:

```
dR/da_i ≈ ( R(u, a_i + h) − R(u, a_i − h) ) / 2h
```

then solves the same linear system with the **same tangent** the hypercomplex path
used, and compares the result with `U^(1)`.

It validates the generated RHS and that solve. It does **NOT** re-solve the
nonlinear problem at perturbed parameters, so it cannot detect an error in the
tangent itself (a wrong `T` appears on both sides and cancels). **It is not a full
solution finite-difference check.** No solution-level FD check is implemented.

## Notes where the config is wider than the code

- `residual.type: element` is accepted, but `resasm_user` handles it with the same
  `PythonResidual` provider as `python` (module + function returning the global
  residual). There is no element-level assembly on this path.
- `sensitivity.backend: dual1` is in the schema, but a **run** on the Python path
  always uses OTILib (`resasm_user/oti_global.py::solve_python`). The value is only
  recorded in the metadata, and `resasm check` rejects `dual1` with `order > 1`.
- `tangent.type: none` is in the schema, but at run time: on the Python path no
  tangent means a hard error; on the black-box path it behaves like `response`.
- The RHS finite-difference check runs on the Python/`element` path only. It is
  never run for `residual.type: executable`.
- On the Python path the real-residual evaluation receives `time` as a float `t`,
  while the OTI order loop receives the pair `(t, t)`. Write your residual to
  tolerate both (most residuals simply ignore `time`).
- On the black-box path `state` is embedded in `request.json`, so it must be
  JSON-serializable there (a `state.file` holding a NumPy array will not serialize).

## Commands

```bash
resasm init --template python|blackbox|cpp|fortran --out DIR   # copy a ready-to-run template
resasm check resasm.yml                                        # readiness report, actionable errors
resasm run   resasm.yml                                        # full sensitivity job
resasm report resasm_output/                                   # summarize a completed run
```

`resasm init --out DIR` refuses to overwrite an existing `DIR` unless you add
`--force`. With no `--template` it falls back to an interactive wizard.

See also: [privacy_contract.md](privacy_contract.md), [glossary.md](glossary.md),
[minimal_user_config.md](minimal_user_config.md).
