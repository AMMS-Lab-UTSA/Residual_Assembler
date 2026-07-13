# Program Definition — Residual_Assembler

> **Residual_Assembler lets a user compute parameter sensitivities from a private
> residual model without sharing their code or mesh. The user only needs to expose
> a local residual evaluator, a converged solution, a tangent, and a small
> configuration file.**

The central object of this program is **not** a UMAT transformation contract.
It is the **residual-sensitivity execution contract**: an agreement about what the
user must be able to evaluate locally, and what the tool will produce from it.

---

## 1. What problem does Residual_Assembler solve?

An engineer has a nonlinear model — an FE code, an in-house solver, a proprietary
material routine — and a converged solution. They want to know **how the solution
responds to the model's parameters**:

- first-order sensitivities `∂u/∂a_i`,
- and, if they ask for it, **higher-order** derivatives `∂²u/∂a_i∂a_j`, `∂³u/…`,
  obtained in the same pass rather than by repeated re-solves.

Finite differences make them re-run the whole nonlinear solve once per parameter
per order, and force them to choose a step size that trades truncation error
against round-off. Hand-differentiating the model is error-prone and has to be
redone whenever the model changes.

The obstacle in practice is usually not the mathematics — it is **secrecy**. The
model, the mesh, and the material law are often the user's intellectual property,
and cannot be handed to a third party or uploaded anywhere.

Residual_Assembler removes both obstacles at once. It computes the sensitivities
from **hypercomplex (OTI/HYPAD) evaluation of the residual**, which yields
arbitrary-order derivatives to machine accuracy in a single pass with no step-size
choice — and it does so **without ever seeing the user's code, mesh, or material
model**. The user keeps the model on their machine and exposes only a residual
evaluator.

## 2. Who is the user?

Someone who already has a working nonlinear model and a converged solution, and
who can evaluate the residual locally. Concretely, one of:

| user | what they have | how they plug in |
|---|---|---|
| a researcher with a Python model | a `residual(u, params)` function | **Path A** — Python provider |
| an engineer with a proprietary solver | a compiled executable they will not share | **Path B** — black-box executable |
| a group with C++/Fortran code | `eval_residual` / `eval_tangent` routines | **Path C** — compiled provider |
| an FE user of this framework | an Abaqus-style model + a registered backend | framework backend (`residual_core/`) |

The user is **not** required to be a developer of this framework, to understand
OTI algebra, or to restructure their model. They are required to be able to
evaluate their own residual.

## 3. What does the user provide?

Five things, and nothing more:

1. **A residual evaluator** — `R(u, a, q, t)`, callable locally.
2. **A converged real solution** `u` — the vector where `R(u, a) ≈ 0`.
3. **A tangent** `T = ∂R/∂u` at that solution — as a function, a file, or returned
   by their executable.
4. **A parameter map** — the names and real values of the parameters of interest.
5. **A small `resasm.yml`** — order, backend, and where the four items above live.

> **The user does not need to provide mesh, element, or material details if their
> residual provider already returns the global residual.** The tool never asks
> what is inside the residual; it only asks that the residual can be evaluated.

Full detail: [input_objects.md](input_objects.md),
[residual_provider_contract.md](residual_provider_contract.md),
[simple_config_contract.md](simple_config_contract.md).

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

The one subtlety that makes higher order work: for `p > 1`, the lower-order
solution sensitivities must be **injected back into the OTI solution vector**
before the next order is evaluated. See
[oti_execution_algorithm.md](oti_execution_algorithm.md) and
[problem_setting.md](problem_setting.md).

## 5. What does the tool output?

A **sensitivity package** in two halves:

- `resasm_output/private/` — the full arrays, for the user only: real residual,
  tangent, `rhs_order<p>.npz`, `solution_sensitivities_order<p>.npz`, dof/parameter
  maps, full validation record.
- `resasm_output/public/` — a report that is **safe to share**: `summary.md`,
  `parameter_ranking.csv`, `sensitivity_norms.csv`, `validation_summary.json`,
  `timing.json`.

Full detail: [output_objects.md](output_objects.md).

## 6. What stays private?

Everything that identifies the model. The `public/` folder never contains source
code, mesh, geometry, the material law, the full residual vector, the tangent
matrix, internal state variables, or private parameter *values* — only parameter
**names** and derived **norms**.

In the black-box path the framework's exposure is structurally limited: it writes
a `request.json`, runs the user's command, and reads a response. It never imports
or inspects the user's code.

Full detail: [privacy_contract.md](privacy_contract.md).

## 7. What is verified?

A ladder of checks, run by `resasm check` (before the job) and `resasm run`
(during it). Levels 0–4 and 6 are implemented today; levels 5, 7 and 8 are defined
but **not** implemented, and are labelled as such.

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

- **We do not solve your nonlinear problem.** You bring a converged solution.
- **We do not verify your solution is converged** unless your residual is exposed.
  In black-box mode the real residual is never returned, so no residual norm can be
  computed; the public report prints `n/a — black-box did not expose real residual`
  and states explicitly that equilibrium was **not** verified. It never prints a
  fabricated `0.000e+00`.
- **We do not check your tangent.** The implemented FD check reuses the same `T` on
  both sides, so an incorrect tangent passes it. A tangent-FD check (level 5) and a
  full solution FD check (level 8) are defined but not implemented.
- **We do not derive your model.** The derivatives are only as correct as the
  residual you expose; if your residual drops the OTI type (e.g. by casting to
  `float`), the imaginary directions are destroyed and the sensitivities are wrong.
  Path A therefore requires generic arithmetic.
- **We add no physics.** `residual_core/` is a separate, formulation-agnostic FE
  residual assembler with explicitly registered backends; it claims support only
  for what is registered and verified. Crystal plasticity is one example backend,
  not the organizing principle.
- **OTILib is not ours and is not vendored.** It is an external GPLv3 library
  (<https://github.com/mauriaristi/otilib.git>) installed separately; on Windows it
  requires WSL.

---

### Where to go next

| you want | read |
|---|---|
| the mathematics | [problem_setting.md](problem_setting.md) |
| the exact algorithm | [oti_execution_algorithm.md](oti_execution_algorithm.md) |
| what to hand the tool | [input_objects.md](input_objects.md) · [residual_provider_contract.md](residual_provider_contract.md) |
| a black-box at order ≥ 2 | [blackbox_order2_contract.md](blackbox_order2_contract.md) — **coefficients, not derivatives** |
| what comes back | [output_objects.md](output_objects.md) |
| the config file | [simple_config_contract.md](simple_config_contract.md) |
| what is checked | [verification_contract.md](verification_contract.md) |
| what stays secret | [privacy_contract.md](privacy_contract.md) |
| the words | [glossary.md](glossary.md) |
| to just run it | [../QUICKSTART_USER.md](../QUICKSTART_USER.md) |
