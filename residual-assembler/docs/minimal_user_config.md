# Minimal user config (`resasm.yml`)

The whole job is described by one small file. Schema:
[`schemas/resasm_config_v1.schema.json`](../schemas/resasm_config_v1.schema.json).

## Python residual (Path A)

```yaml
problem:
  name: spring_demo
  unknowns: 1

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

## Black-box executable (Path C)

```yaml
problem:
  name: private_solver_demo
  unknowns: 1500

residual:
  type: executable
  command: ./my_solver --request {request} --response {response}

tangent:
  type: file
  file: tangent.npz

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

## Fields

| Field | Required | Meaning |
|---|---|---|
| `problem.name` | yes | a label for the run |
| `problem.unknowns` | **no** | number of DOFs. Omit it and it is inferred from the solution vector. If you do give it, it must match the solution length (a mismatch is an error, never a silent truncation). |
| `residual.type` | yes | `python` \| `executable` \| `element` |
| `residual.module` | python/element | path to your `.py` |
| `residual.function` | no | residual function name (default `residual`) |
| `residual.command` | executable | command with `{request}` / `{response}` |
| `tangent.type` | no | `python` \| `file` \| `response` \| `none` |
| `tangent.function` / `tangent.file` | as needed | tangent source |
| `parameters` | yes | `name: value` design parameters (seeded) |
| `solution.file` or `solution.values` | yes | converged `u` |
| `sensitivity.order` | yes | derivative order (≥ 1) |
| `sensitivity.backend` | no | `otilib` (default) or `dual1` (order 1) |
| `constraints.free` / `constraints.prescribed` | no | free/fixed DOF indices (default: all free) |
| `state.file` / `time` | no | history/time passed to your residual |
| `output.dir` | no | output folder (default `resasm_output`) |
| `validation.rhs_finite_difference_check` | no | finite-difference the **residual derivative** `d(residual)/d(parameter)` at the fixed `u` and solve with the same tangent, then compare to the order-1 sensitivity. Validates the generated RHS + that solve. It does **not** re-solve the nonlinear problem, so it cannot catch an error in the tangent itself. (Python path.) |
| `validation.finite_difference` | no | deprecated alias of the above — same behaviour, misleading name. |
| `validation.solution_finite_difference_solver` | no | **reserved, not implemented.** A true solution-level FD check would re-run your nonlinear solver at perturbed parameters. Setting it performs no such check; the run says so explicitly. |

No mesh/element/material fields are required when your residual returns the
global residual. See [user_input_contract.md](user_input_contract.md).
