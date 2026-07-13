# Output Contract

The verification zoo is a proving ground; **the deliverable is a residual /
sensitivity package** that supplies the numerical ingredients of the HYPAD
residual method (Aristizabal et al., *Non-Intrusive Arbitrary-order Sensitivity
Analysis of Non-Linear Residual-Based Problems*):

$$\mathbf{T}\,\mathbf{U}^{(p)} = -\mathbf{R}^{(p)}$$

where $\mathbf{T}$ is the tangent at the converged real solution, $\mathbf{R}^{(p)}$
is assembled from the $p$-th order imaginary coefficients of a hypercomplex
residual evaluation, and $\mathbf{U}^{(p)}$ are the unknown $p$-th order solution
derivative coefficients. The framework answers three questions:

1. **Can I assemble this residual?**
2. **What minimum data is missing?**
3. **If yes, what are $\mathbf{R}$, $\mathbf{T}$ and the RHS matrices for sensitivities?**

> **OTILib is the production HYPAD backend.** Arbitrary-order, multi-parameter
> $\mathbf{R}^{(p)}$ is produced by the OTILib provider
> (`core/oti_rhs_provider.py` + `algebra/otilib_adapter.py`), which seeds all
> parameters simultaneously ($a_i^* = a_i + e_i$) and solves order-by-order,
> injecting the solved lower-order $\mathbf{U}^{(k)}$ into $\mathbf{u}^*$ before
> advancing to order $p$. The first-order `Dual1` provider was only a preliminary
> smoke test. The contract below fixes the shapes, orderings and maps so any
> hypercomplex residual evaluation drops its $\mathbf{R}^{(p)}$ columns in — see
> §4 and `set_residual_order`. When OTILib is not installed the OTI backend is
> reported cleanly (no silent fallback).

## Objects (in `residual_core.core`)

| Object | Module | Role |
|---|---|---|
| `ResidualAssemblyResult` | `results.py` | real residual $\mathbf{R}$, reactions, element contributions, embedded tangent |
| `TangentResult` | `results.py` | tangent $\mathbf{T}$ + source + DOF partition |
| `StateResult` | `results.py` | state variables / sensitivities / history status |
| `ValidationReport` | `results.py` | residual/reaction/UMAT/UEL comparisons + tolerances |
| `SensitivityRHSResult` | `sensitivity_package.py` | $\mathbf{R}^{(p)}$, $\mathbf{rhs}^{(p)}=-\mathbf{R}^{(p)}$, maps, algebra metadata |
| `AlgebraMetadata` | `sensitivity_package.py` | OTI basis count $m$, truncation order $n_t$, counts |
| `SensitivityPackage` | `sensitivity_package.py` | the whole deliverable; `solve(p)` |

Public entry points: `ResidualProblem.result(mode=...)` and
`ResidualProblem.sensitivity_package(mode=..., parameters=..., max_order=...)`.

## Serialization

Every result implements a common contract:

- `to_dict()` → JSON-serializable **metadata** (norms, shapes, maps, flags),
- `arrays()` → named `numpy` arrays,
- `to_markdown()` → human-readable report,
- `save(prefix)` → writes `prefix.json` + `prefix.npz` + `prefix.md`.

`SensitivityPackage.save(prefix)` additionally writes `prefix_residual.*`,
`prefix_sensitivity.*`, `prefix_validation.*`, `prefix_state.*`.

---

## 1. Inspection output

From `ResidualProblem.inspect()` / embedded in the package as `inspection`:

- **supported backends** — per element type, the selected formulation backend and
  its verification status;
- **available modes** — which assembly modes are reachable for the model;
- **minimum missing data** — the single next input required, never a generic
  checklist.

## 2. Real residual output — `ResidualAssemblyResult`

- **global residual vector** `R` (array, length `ndof`);
- **free-DOF residual norm** `||R_free||` (≈ 0 at a converged solution);
- **constrained-DOF reaction comparison** when an exported `RF` is attached
  (`reaction_comparison()`: relative error, detected sign, pass/fail);
- **element residual contributions** — captured as a CSR-like triple
  (`element_ids`, `element_res_ptr`, `element_res_dofs`, `element_res_values`);
- **sign convention** — `R = F_internal − F_external`; Abaqus reports `RHS = −R`.

## 3. Tangent output — `TangentResult`

- **tangent matrix** `T` (array) when a backend assembles one;
- **source** — `solver-exported` | `backend-assembled` | `finite-difference` |
  `hypad` | `unavailable`;
- **DOF ordering** — global labels `node:TYPE` (e.g. `12:RZ`) in assembly order;
- **free / prescribed partition** — `free_mask`, `prescribed_idx`.

A mode with no material/element tangent (e.g. stress-driven) reports
`source = unavailable` and an all-zero `T` is **never** presented as valid — the
package's `tangent_available()` is `False` and `solve()` refuses with guidance.

## 4. Sensitivity RHS output — `SensitivityRHSResult`

- **raw** `R^(p)` — shape `(ndof, N^(p))`, `N^(p) = C(p+m−1, p)`; one column per
  order-$p$ imaginary direction;
- **`rhs^(p) = -R^(p)`** — the right-hand-side of the linear system;
- **derivative-direction map** — per order, each column's exponent multi-index,
  human label (`e1^2`, `e1*e2`), the participating parameters, and the
  **recovery factor** $\prod \kappa_i!$ to convert the extracted imaginary
  coefficient into the true partial derivative;
- **parameter map** — design-parameter name → imaginary basis index (`e_i`);
- **order `p`** and **algebra metadata** — algebra name, basis count $m$,
  truncation order $n_t$, total coefficient count $N = C(m+n_t, m)$.

Plug-in point: `set_residual_order(p, array)` accepts an externally-computed
$\mathbf{R}^{(p)}$ (from an OTI/HYPAD residual evaluation) with shape-checking; the
rest (rhs, solve, export) then works unchanged.

### Who generates $\mathbf{R}^{(p)}$? (stress-driven vs hypercomplex)

There are two distinct capabilities, and they must not be conflated:

- **stress-driven mode gives $\mathbf{R}$** — assembling the residual from an
  *exported* field verifies the finite-element assembly, but it **cannot** produce
  a parameter-sensitivity right-hand-side by itself. The exported field is a fixed
  set of numbers with no dependence on the design parameters, so $\partial R /
  \partial a_i$ is not recoverable from it (unless *field sensitivities* are also
  exported).
- **material/formulation hypercomplex mode gives $\mathbf{R}^{(p)}$** — a
  *parameterized* backend whose residual can be evaluated hypercomplexly (with a
  perturbed parameter) yields $\partial R/\partial a_i$ directly from the imaginary
  part. This is what an **RHS provider** does.

An **RHS provider** (`core/rhs_provider.py`) is what actually fills
$\mathbf{R}^{(p)}$:

```python
SensitivityRHSProvider.evaluate_rhs(problem, real_solution, tangent,
                                    parameters, order) -> SensitivityRHSResult
```

`DualNumberRHSProvider` (order 1, **legacy**) perturbs one parameter at a time
along a dual imaginary direction, evaluates each element residual at the converged
real solution, and reads $\mathbf{R}^{(1)}[:,i] = \mathrm{Im}_i[r(u,
a+\varepsilon_i)] = \partial R/\partial a_i$. It reports `hypercomplex_ready =
False` (leaving the column zero) for backends whose residual is not dual-safe — the
contract degrades cleanly instead of fabricating a derivative. The first verified
provider + backend pair is `DualNumberRHSProvider` + `nonlinear_spring1`
(`tests/framework/test_nonlinear_spring_sensitivity.py`).

`OtiLibRHSProvider` (`core/oti_rhs_provider.py`) is the **production** provider:
one OTILib evaluation carries all $m$ parameters (seeded $a_i^* = a_i + e_i$) to
truncation order $n_t = q$. It loops $p = 1 \dots q$, extracts the order-$p$
imaginary coefficients into $\mathbf{R}^{(p)}$, solves $\mathbf{T}\,\mathbf{U}^{(p)}
= -\mathbf{R}^{(p)}$ on the free partition, and **injects** $\mathbf{U}^{(p)}$ back
into $\mathbf{u}^*$ before advancing — the order-by-order coupling required for
$p > 1$. It requires a genuine OTILib install and reports cleanly when absent
(no fallback to Dual1). CLI: `resasm sensitivity model.json --params params.json
--order q --backend otilib`.

Use it via the facade:

```python
pkg = problem.sensitivity_package(mode="formulation", parameters=["spring.k"],
                                  generate_rhs=True, fd_check=True)
U1  = pkg.solve(1)     # real du/dk, cross-checked against finite differences
```

or the CLI: `resasm sensitivity model.json --params params.json --order 1`.

## 5. State / history output — `StateResult`

- **state variables** (per element / integration point);
- **state sensitivities** when available;
- **load increment index**, **time**, **dtime**;
- **history replay status** — for history-dependent materials, whether the full
  increment sequence is present (final-step-only replay is invalid).

## 6. Validation output — `ValidationReport`

A table of checks, each with quantity, value, tolerance and pass/fail:

- **residual-vs-solver** — `||R_free||` at the solver solution;
- **reaction-force** — computed reactions vs exported `RF` (sign-aware);
- **UMAT replay** — `STRESS` / `STATEV` vs ODB (history in order);
- **UEL RHS/AMATRX** — `element_residual = −RHS`, `AMATRX = dR/du`.

Pending checks (data not yet available) are reported as `pending`, never as pass.

---

## The system, assembled

```python
from residual_core import ResidualProblem

p   = ResidualProblem.from_abaqus("model.inp")
pkg = p.sensitivity_package(mode="material-replay",
                            parameters=["mat.E", "mat.nu"], max_order=3)

pkg.runnable            # can I assemble?  (else pkg.minimum_missing)
pkg.residual.R          # R (order 0)
pkg.residual.tangent.T  # T
pkg.sensitivity.expected_shape(2)      # shape of R^(2)
pkg.sensitivity.direction_map(2)       # columns ↔ derivative directions

# OTILib backend fills the RHS order-by-order; or drop columns in directly:
pkg.sensitivity.set_residual_order(2, R2)   # R2 from an OTILib residual evaluation
U2 = pkg.solve(2)                            # T U^(2) = -R^(2)

pkg.save("out/model_sensitivity")            # .json + .npz + .md bundle
```

Prescribed DOFs (boundary conditions independent of the design parameters) carry
zero sensitivity, so `solve` works on the free partition and scatters zeros back.
A properly-constrained, converged FE problem yields a non-singular free tangent;
an unconstrained model (rigid-body modes, or a pin-jointed truss) is rejected with
a clear message rather than a silent wrong answer.
