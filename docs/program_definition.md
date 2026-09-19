# Program definition: Residual_Assembler

This page defines the program precisely: the problem it solves, who uses it,
what they provide, what it computes and returns, what stays private, what is
verified and what is not claimed. It is for readers who need the exact scope;
a one-page summary is [definition_outline.md](definition_outline.md).

> **Residual_Assembler computes parameter sensitivities of a converged
> nonlinear model without re-solving it and without the user sharing code,
> mesh or material law.** For a finished Abaqus analysis it needs the input
> deck, the output database, the material compiled as an OTI provider and a
> request. For a model whose residual the user can evaluate, it needs a local
> residual evaluator, a converged solution, a tangent and a small
> configuration file.

The central object of the program is the **residual-sensitivity execution
contract**: an agreement about what the user must make available locally and
what the tool will produce from it. For an analysis replay the contract is
the four files and the provider's `Mapping.json`
([REQUEST_INTERFACE.md](REQUEST_INTERFACE.md)); for a supplied residual it is
`resasm.yml` and the provider interface described below.

---

## 1. What problem does Residual_Assembler solve?

An engineer has a nonlinear model (a finite-element analysis, an in-house
solver, a proprietary material routine) and a converged solution. They want to
know **how the solution responds to the model's parameters**:

- first-order sensitivities `∂u/∂a_i`, and of derived results such as
  reactions, stresses and state variables;
- and, where the residual is supplied directly, **higher-order** derivatives
  `∂²u/∂a_i∂a_j`, `∂³u/…`, obtained in the same pass rather than by repeated
  re-solves.

Finite differences make them re-run the whole nonlinear solve once per
parameter per order, and force them to choose a step size that trades
truncation error against round-off. Hand-differentiating the model is
error-prone and has to be redone whenever the model changes.

The obstacle in practice is usually not the mathematics but **secrecy**. The
model, the mesh and the material law are often the user's intellectual
property, and cannot be handed to a third party or uploaded anywhere.

Residual_Assembler removes both obstacles. It computes the sensitivities from
a **hypercomplex (OTI/HYPAD) evaluation of the residual**, which yields
derivatives to machine accuracy with no step size to choose, and it does so
**without seeing the user's code**. For an Abaqus analysis the material
routine arrives compiled, as an OTI provider built by the companion UMAT-OTI;
for a supplied residual the user exposes only an evaluator.

## 2. Who is the user?

Someone who already has a working nonlinear model and a converged solution.
Concretely, one of:

| user | what they have | how they plug in |
|---|---|---|
| an Abaqus user with a finished analysis | `Analysis.inp`, `Analysis.odb`, and a UMAT built into an OTI provider | **analysis replay**: `resasm request` / `resasm history` |
| an FE user who wants the residual itself | a mesh, exported fields, a material | **Path A**: assembly from ingredients |
| an engineer with a proprietary solver | a compiled executable they will not share | **Path B**: black-box executable |
| a group with C++/Fortran code | `eval_residual` / `eval_tangent` routines | **Path B**: compiled black-box provider |
| a researcher with a Python model | a `residual(u, params)` function | **Path C**: direct residual |

The user is **not** required to be a developer of this framework, to
understand OTI algebra, or to restructure their model.

## 3. What does the user provide?

For an **analysis replay**: the input deck, the output database, the OTI
provider object with its `Mapping.json`, and a `sensitivity_request.json`
naming the outputs, parameters, region and increments. The provider is built
once by whoever owns the UMAT source; the person running the replay does not
need that source.

For a **supplied residual**, five things, and nothing more:

1. **A residual evaluator** — `R(u, a, q, t)`, callable locally.
2. **A converged real solution** `u` — the vector where `R(u, a) ≈ 0`.
3. **A tangent** `T = ∂R/∂u` at that solution — as a function, a file, or returned
   by their executable.
4. **A parameter map** — the names and real values of the parameters of interest.
5. **A small `resasm.yml`** — order, backend, and where the four items above live.

> **The user does not need to provide mesh, element, or material details if their
> residual provider already returns the global residual.** The tool never asks
> what is inside the residual; it only asks that the residual can be evaluated.

Full detail: [REQUEST_INTERFACE.md](REQUEST_INTERFACE.md) (analysis replay);
[input_objects.md](input_objects.md),
[residual_provider_contract.md](residual_provider_contract.md),
[simple_config_contract.md](simple_config_contract.md) (supplied residual).

## 4. What does the tool compute?

At the converged solution the residual is stationary in the solution variable, so
differentiating `R(u(a), a) = 0` with respect to the parameters gives the
**sensitivity equation** solved order by order:

```
    T U^(p) = -R^(p)
```

| symbol | meaning | where it comes from |
|---|---|---|
| `T`     | tangent/Jacobian `∂R/∂u` at the converged real solution | the user's tangent provider |
| `R^(p)` | p-th order residual coefficient matrix | OTI evaluation of the user's residual |
| `U^(p)` | p-th order solution sensitivity coefficients | the linear solve above |

The tool does **not** solve the user's nonlinear problem. The real solve is done
first, by the user. Residual_Assembler re-evaluates the *already converged*
residual with OTI-perturbed parameters and extracts the derivative coefficients.

For a history-dependent analysis the equation is solved at every recorded
increment, `K du_n/dp = -dR_n/dp`, with the stress and state derivatives
carried from one increment to the next; see
[REPLAY_HISTORY.md](REPLAY_HISTORY.md).

The one subtlety that makes higher order work: for `p > 1`, the lower-order
solution sensitivities must be **injected back into the OTI solution vector**
before the next order is evaluated. See
[oti_execution_algorithm.md](oti_execution_algorithm.md) and
[problem_setting.md](problem_setting.md).

## 5. What does the tool output?

An **analysis replay** writes `sensitivity_results.json`,
`sensitivity_tables.csv` and `run_report.txt` at the top of its output
directory, and keeps full fields, `K`, `R` and derivatives in `private/`.

A **supplied-residual run** writes a sensitivity package in two halves:

- `resasm_output/private/` — the full arrays, for the user only: real residual,
  tangent, `rhs_order<p>.npz`, `solution_sensitivities_order<p>.npz`, dof/parameter
  maps, full validation record.
- `resasm_output/public/` — a report that is **safe to share**: `summary.md`,
  `parameter_ranking.csv`, `sensitivity_norms.csv`, `validation_summary.json`,
  `timing.json`.

Full detail: [REQUEST_INTERFACE.md](REQUEST_INTERFACE.md),
[output_objects.md](output_objects.md).

## 6. What stays private?

Everything that identifies the model. No command makes a network call. The
`public/` folder of a supplied-residual run never contains source code, mesh,
geometry, the material law, the full residual vector, the tangent matrix,
internal state variables, or private parameter *values*, only parameter
**names** and derived **norms**. An analysis replay keeps every full array in
`private/`, and needs no material source: the provider arrives compiled.

In the black-box path the framework's exposure is structurally limited: it writes
a `request.json`, runs the user's command, and reads a response. It never imports
or inspects the user's code.

Full detail: [privacy_contract.md](privacy_contract.md).

## 7. What is verified?

An **analysis replay** checks the replayed stress, state and reactions against
the ODB at every integration point and increment, and fails rather than warns
on any excess. On request it also checks the provider tangent and compares the
sensitivities with whole-model central differences of the ORIGINAL UMAT
(`resasm history --verify tangent|fd`, `resasm request --validate`). The
measured results are in [VERIFICATION_RECORD.md](VERIFICATION_RECORD.md).

For a **supplied residual**, a ladder of checks runs in `resasm check` (before
the job) and `resasm run` (during it). Levels 0–4 and 6 are implemented today;
levels 5, 7 and 8 are defined but **not** implemented, and are labelled as
such.

The one that matters most for honesty: the finite-difference check that *is*
implemented is a **RHS finite-difference check** — it differences
`d(residual)/d(parameter)` at the **fixed** solution `u` and solves with the
**same** tangent. It validates the generated right-hand side and the solve.
**It is not a full solution finite-difference validation**, it does not re-solve
the nonlinear problem, and it cannot detect an error in the tangent itself
(a wrong `T` cancels from both sides).

Full detail: [verification_contract.md](verification_contract.md).

## 8. What is not claimed?

Stated plainly, so nobody infers more than is true:

- **We do not solve your nonlinear problem.** You bring a converged solution or
  a converged analysis. (`resasm history --reequilibrate` Newton-polishes each
  recorded increment to double-precision equilibrium, starting from the
  recorded state; it does not re-run the analysis.)
- **We do not verify your solution is converged** unless your residual is exposed.
  In black-box mode the real residual is never returned, so no residual norm can be
  computed; the public report prints `n/a — black-box did not expose real residual`
  and states explicitly that equilibrium was **not** verified. It never prints a
  fabricated `0.000e+00`.
- **We do not check a tangent you supply.** On the supplied-residual path the
  implemented FD check reuses the same `T` on both sides, so an incorrect
  tangent passes it. A tangent-FD check (level 5) and a full solution FD check
  (level 8) are defined but not implemented there. (The analysis replay builds
  its own tangent from the provider and can check it: `--verify tangent`.)
- **We do not derive your model.** The derivatives are only as correct as the
  residual you expose; if your residual drops the OTI type (e.g. by casting to
  `float`), the imaginary directions are destroyed and the sensitivities are wrong.
  A Python residual (Path C) therefore requires generic arithmetic.
- **We add no physics.** `residual_core/` is a formulation-agnostic FE residual
  assembler with explicitly registered backends; it claims support only for
  what is registered and verified. The analysis replay supports small-strain
  C3D8 analyses with one static step and one user material, and refuses
  anything else by name ([STATUS.md](../STATUS.md)). Crystal plasticity is one
  example backend, not the organizing principle.
- **OTILib is not ours and is not vendored.** It is an external GPLv3 library
  (<https://github.com/mauriaristi/otilib.git>) installed separately; on Windows it
  requires WSL.

---

### Where to go next

| you want | read |
|---|---|
| to compute sensitivities of an Abaqus analysis | [REQUEST_INTERFACE.md](REQUEST_INTERFACE.md) · [REPLAY_HISTORY.md](REPLAY_HISTORY.md) |
| the mathematics | [problem_setting.md](problem_setting.md) |
| the exact algorithm | [oti_execution_algorithm.md](oti_execution_algorithm.md) |
| what to hand the tool | [input_objects.md](input_objects.md) · [residual_provider_contract.md](residual_provider_contract.md) |
| a black-box at order ≥ 2 | [blackbox_order2_contract.md](blackbox_order2_contract.md) — **coefficients, not derivatives** |
| what comes back | [output_objects.md](output_objects.md) |
| the config file | [simple_config_contract.md](simple_config_contract.md) |
| what is checked | [verification_contract.md](verification_contract.md) |
| what stays secret | [privacy_contract.md](privacy_contract.md) |
| the words | [glossary.md](glossary.md) |
| to just run it | [../QUICKSTART_USER.md](../QUICKSTART_USER.md) · [CLI_GUIDE.md](CLI_GUIDE.md) · [../examples/README.md](../examples/README.md) |
