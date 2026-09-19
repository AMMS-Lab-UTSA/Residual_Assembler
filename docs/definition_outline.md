# Residual_Assembler in twelve points

A condensed definition of the program: the problem it solves, the equation at
its centre, what users provide, what it computes and returns, how privacy and
verification work, and where it stands. It is for readers who want the whole
picture on one page before reading the detailed contracts; every point links
to the page that specifies it. The longer statement is
[program_definition.md](program_definition.md).

## 1. The problem

- A user has a converged nonlinear model, often a finite-element analysis, and
  wants the parameter sensitivities `du/da` of its results.
- The model, mesh and material are often private: they cannot be sent out or
  read by an outside tool.
- Finite differences need at least one extra nonlinear solve per parameter,
  and the step size is a guess.
- Residual_Assembler differentiates the **residual** at the converged solution
  instead. It only needs to evaluate that residual, never to read the code
  behind it, and it obtains exact first- and higher-order sensitivities without
  re-solving the model.

## 2. The residual sensitivity equation

At the converged real solution `R(u, a) = 0`. Differentiating with respect to
the parameters gives one linear system per order `p`:

```
T U^(p) = -R^(p)
```

- `T` is the tangent `dR/du` at the converged solution.
- `R^(p)` is the order-`p` residual coefficient matrix from the OTI
  (order-truncated imaginary) evaluation of the residual.
- `U^(p)` is the order-`p` solution sensitivity: the answer.
- The same `T` serves every order: factor it once and change only the
  right-hand side. The system is solved on the free DOFs; prescribed DOFs have
  zero sensitivity.

For a history-dependent finite-element analysis the same equation is solved at
every increment, `K du_n/dp = -dR_n/dp`, with the derivatives of stress and
state carried from one increment to the next (point 6).

Details: [problem_setting.md](problem_setting.md).

## 3. Three ways to obtain the residual

| Route | The user provides | Command |
|---|---|---|
| **Replay a finished analysis** | `Analysis.inp`, `Analysis.odb`, the compiled OTI material provider `OTI_UMAT.obj` + `Mapping.json`, `sensitivity_request.json` | `resasm request`, `resasm history` |
| **Assemble from ingredients** (Path A) | mesh, element formulation, material, solution field, loads and constraints | `resasm inspect`, `init-assembly`, `assemble`, `run` |
| **Supply the residual** | a Python `residual(u, params)` (Path C) or a private executable that returns residual coefficients (Path B) | `resasm init`, `check`, `run` |

Abaqus never exposes its global residual, so an Abaqus user takes the first
or second route: Residual_Assembler rebuilds `R` from what the analysis does
export. Guidance: [which_path_should_i_use.md](which_path_should_i_use.md).

## 4. What the user provides

For the analysis replay: the four files above. The material developer builds
the provider once with the companion UMAT-OTI and shares the compiled object
and its `Mapping.json`; the UMAT source is not needed to run the replay.

For a supplied residual (`resasm.yml`):

- a converged solution vector `u` (`solution.file` or `solution.values`);
- parameter names and real values (`parameters: {k: 2.0, f: 16.0}`), the design
  variables that are seeded;
- a residual evaluator: a Python function, an executable, or compiled C++ or
  Fortran behind the same request/response contract;
- a tangent `T`: a Python `tangent()`, a saved `tangent.npz`, or returned in the
  black-box response;
- the derivative order, `sensitivity.order`.

No mesh, element data or material model is needed when the evaluator returns
the global residual. Details: [input_objects.md](input_objects.md),
[residual_provider_contract.md](residual_provider_contract.md).

## 5. The OTI execution loop

For a supplied residual (`resasm_user/oti_global.py::solve_python`):

1. **Define the number system:** `m` parameters, truncation order `q`. One OTI
   number carries every parameter direction at once.
2. **Seed** every parameter at once: `a_i* = a_i + e_i`.
3. **Evaluate** `R* = residual(u*, a*)`, the user's function called with OTI
   scalars.
4. **Extract** `R^(p)`, the order-`p` coefficients of `R*`, one column per
   direction (including mixed directions such as `e1*e2`).
5. **Solve and re-inject:** `U^(p) = -T_ff^-1 R^(p)` on the free DOFs, then
   write `U^(p)` back into `u*` **before** order `p+1`.

The loop runs `p = 1 … q`: one residual evaluation per order, not one per
parameter. The direction list of each order and each direction's recovery
factor (which converts an OTI coefficient into the true partial derivative)
come from `residual_core/algebra/otilib_adapter.py::order_directions`.
Details: [oti_execution_algorithm.md](oti_execution_algorithm.md).

## 6. History-dependent analyses

For plasticity, crystal plasticity and other path-dependent materials the
sensitivity at an increment depends on the whole history before it, so the
replay marches every recorded increment. At each increment and integration
point, one call of the provider's `UMAT_OTI_EVAL_TOTAL` returns the stress,
state and tangent with their parameter derivatives, seeded with the derivatives
carried from the previous increment. The engine assembles `K` and `dR/dp`
sparsely, factorises `K_ff` once per increment for all parameters, and solves
for `du/dp`; reactions, stresses, von Mises stress and state follow. The
replayed stress, state and reactions are checked against the ODB at every
integration point and increment. Details: [REPLAY_HISTORY.md](REPLAY_HISTORY.md).

## 7. The output package

- **Analysis replay:** `sensitivity_results.json`, `sensitivity_tables.csv` and
  `run_report.txt` at the top of the output directory; full fields, `K`, `R`
  and derivatives in `private/` ([REQUEST_INTERFACE.md](REQUEST_INTERFACE.md)).
- **Supplied residual:** one directory per run (default `resasm_output/`),
  split into `private/` (`metadata.json`, `parameter_map.json`, `dof_map.json`,
  `residual_real.npz`, `tangent.npz`, `rhs_order<p>.npz`,
  `solution_sensitivities_order<p>.npz`, `direction_map_order<p>.json`,
  `validation_full.json`) and `public/` (`summary.md`, `parameter_ranking.csv`,
  `sensitivity_norms.csv`, `validation_summary.json`, `timing.json`).
  `resasm report resasm_output/` prints the headline; `from resasm_user import
  read_report` reads it in Python ([output_objects.md](output_objects.md)).

## 8. Privacy model

- The tool calls the residual or the provider; it never reads its source. On
  the black-box path it only exchanges two files with an executable it knows
  nothing about. Nothing is sent over the network.
- **Never required:** source code; on the direct and black-box paths also mesh,
  element connectivity, material models and state variables.
- **Never in `public/`:** source, mesh, the full residual vector, the full
  tangent, state variables, or parameter *values*; only parameter **names**
  and sensitivity **norms**. This is enforced where the outputs are written.
- **`private/`** holds every full array and stays on the user's disk unless
  they choose to send it.
- On the black-box path the real residual is never exposed, so
  `residual_free_norm` is written as `null` and `summary.md` says "Equilibrium
  was NOT verified by this run", never a fabricated `0.0`.

Details: [privacy_contract.md](privacy_contract.md).

## 9. Verification

- **Analysis replay:** stress, state and reactions against the ODB at every
  integration point and increment; optionally the provider tangent and
  whole-model central differences of the ORIGINAL UMAT
  (`resasm history --verify tangent|fd`, `resasm request --validate`); and, for
  suitable models, Euler's homogeneity identity at every increment. Every
  quantitative claim, with its command and measured value, is in
  [VERIFICATION_RECORD.md](VERIFICATION_RECORD.md).
- **Supplied residual:** `resasm check resasm.yml` is a step-by-step readiness
  report (`[ok]`, `[warn]` or `[fail]` per line) that stops at the first
  actionable blocker and shows the YAML to add. The RHS finite-difference
  cross-check (`validation.rhs_finite_difference_check`) central-differences
  `dR/da_i` at the fixed `u`, solves with the same tangent and compares to
  `U^(1)`; it passes below a relative error of 1e-4. It does **not** re-solve
  the nonlinear problem, so a wrong tangent cancels on both sides, and the
  report says so. `validation.solution_finite_difference_solver` is reserved
  and not implemented.
- **Assembly framework:** a generic ladder of Levels 0 to 7 (zero field, rigid
  body, patch test, FD tangent, solver comparison, material replay, full
  residual replay) in `residual_core/core/verification.py`.

Details: [verification_contract.md](verification_contract.md).

## 10. Example: a one-line spring residual

`templates/user_python_residual/` holds the whole model in two functions:

```python
def residual(u, params, state=None, time=None):
    return [params["k"] * u[0] ** 3 - params["f"]]     # R = k u^3 - f

def tangent(u, params, state=None, time=None):
    return [[3.0 * params["k"] * u[0] ** 2]]
```

With `parameters: {k: 2.0, f: 16.0}`, `solution.npy = [2.0]` and
`sensitivity: {order: 2, backend: otilib}`, the run reproduces the analytic
`du/dk = -u/(3k) = -1/3` and `du/df = 1/(3k u^2) = 1/24`, and the RHS
cross-check confirms order 1. Order 2 comes from the same run, because the OTI
number already carries `e1^2`, `e1*e2` and `e2^2`; the recovered
`d²u/dk² = 2/9`. Tested in `tests/framework/test_otilib_spring_sensitivity.py`.

## 11. Example: a private black-box solver

`templates/user_blackbox_residual/` keeps the model inside the user's
executable. `resasm.yml` says only
`residual: {type: executable, command: python my_solver.py --request {request} --response {response}}`
and `tangent: {type: response}`. For each order `p` the framework writes
`request.json` (`u`, `parameters`, `direction_map`, `u_star_coefficients`, …)
and reads back `response.npz` with `R_order_<p>` and, optionally, `tangent`.
The executable differentiates; the framework solves. For order 2 and above the
executable must use `u_star_coefficients` (the lower-order `U^(k)` already
solved) to rebuild `u*`, and must return Taylor coefficients, not derivatives
([blackbox_order2_contract.md](blackbox_order2_contract.md)). OTILib is not
needed on the framework side for this path.

## 12. Status and next steps

Supported and verified today ([STATUS.md](../STATUS.md)):

- total-history sensitivities of small-strain C3D8 Abaqus analyses for any
  UMAT-OTI provider, including two full-size cantilevers (J2, 1,536 C3D8 and
  40 increments; FCC crystal plasticity, 384 C3D8 and 25 increments), checked
  against whole-model finite differences of the ORIGINAL UMAT, Abaqus reruns
  and the homogeneity identity;
- stress-driven C3D8 assembly, and bounded finite-strain neo-Hookean assembly
  with first-order parameter sensitivities;
- the direct Python and black-box residual paths at arbitrary order, with five
  templates (`python`, `blackbox`, `blackbox-order2`, `cpp`, `fortran`);
- a clean-install gate that reproduces the workflows from wheels built from
  clean clones.

OTILib is external (GPLv3), not vendored, and does not build natively on
Windows; there the Python residual path needs WSL. The OTILib tests pass with
a genuine build (11 passed, none skipped, with `RUN_OTILIB_TESTS=1`).

Open items:

- finite-strain plasticity, element types other than C3D8, several steps or
  materials, and distributed loads;
- the Oxford crystal-plasticity UMAT, which needs Intel ifort to compile;
- the finite-strain `DDSDDE` to `AMATRX` mapping of `solid_c3d8_finite_strain`;
- OTI differentiation through the assembly recipe for C3D8 backends (today
  C3D8 sensitivities come from the replay engines and the bounded
  finite-strain path);
- the reserved solution-level finite-difference check for supplied residuals,
  which would re-run the user's solver at perturbed parameters, the one check
  that catches a wrong user tangent.
