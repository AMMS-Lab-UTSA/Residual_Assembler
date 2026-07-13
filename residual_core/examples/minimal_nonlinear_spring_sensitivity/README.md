# Minimal example — nonlinear spring parameter sensitivity

The first **real** vertical slice of the residual-sensitivity method: a backend
that actually *generates* the right-hand-side $\mathbf{R}^{(1)}$ from a
hypercomplex (dual-number) residual evaluation, then solves

$$\mathbf{T}\,\mathbf{U}^{(1)} = -\mathbf{R}^{(1)}$$

for a parameter sensitivity — verified against the closed-form answer and finite
differences.

## The toy problem

A single grounded degree of freedom with a cubic force law
(`formulations/nonlinear_spring1.py`):

$$R(u,k,f) = k\,u^3 - f$$

At equilibrium $u = (f/k)^{1/3}$, with

$$T = \frac{\partial R}{\partial u} = 3k\,u^2,\qquad
\frac{\partial R}{\partial k} = u^3,\qquad
\frac{du}{dk} = \frac{-\partial R/\partial k}{T} = -\frac{u}{3k}.$$

Model: [model.json](model.json) — `k = 2`, `f = 16` ⟹ $u = (16/2)^{1/3} = 2$, so
$du/dk = -2/(3\cdot2) = -1/3$.

## Command

```bash
resasm sensitivity residual_core/examples/minimal_nonlinear_spring_sensitivity/model.json \
       --params residual_core/examples/minimal_nonlinear_spring_sensitivity/params.json --order 1
```

(Use `python -m residual_core.ui.cli sensitivity ...` if the console script is not
installed.)

## Expected output

```
sensitivity analysis (mode=formulation, order=1)
  residual norm ||R_free|| = 0.000000e+00
  tangent source           = backend-assembled
  parameters               = ['spring.k']
  R^(1) shape              = (1, 1)
  hypercomplex ready       = True
  solved sensitivities du/da:
    spring.k       : -3.333333e-01   [analytic -3.333333e-01, rel 0.00e+00]   [FD -3.333330e-01, rel 1.00e-06]
```

- **residual norm** ≈ 0 confirms the Newton solve converged (`R(u) = 0`).
- **R^(1) = u³ = 8** was produced by evaluating the residual with a *dual* value
  of `k` and reading `dR/dk` from the imaginary part — no hand-coded derivative.
- **du/dk = −1/3** matches the analytical value exactly and the finite-difference
  value to ~1e-6.

## What was auto-detected vs supplied

- Auto-detected: element type `SPRING1` → `nonlinear_spring1` backend; the
  assembly mode (`formulation`); the tangent (backend-assembled); the design
  parameter list is read from the material section keys if not given.
- Supplied: the parameter of interest (`spring.k`) and the section values
  (`k`, `f`) in the model; optionally an `expected` analytic value in
  `params.json` for the printed comparison.

## Python API

```python
from residual_core import ResidualProblem
p   = ResidualProblem.from_neutral(".../model.json")
p.solve_newton("formulation")                       # real solution u
pkg = p.sensitivity_package(mode="formulation", parameters=["spring.k"],
                            generate_rhs=True, fd_check=True)
U1  = pkg.solve(1)                                  # du/dk = -u/(3k)
pkg.save("out/spring_sensitivity")                  # .json + .npz + .md
```

## Important limitation

Stress-driven residual assembly (from an exported field) yields `R` but **cannot**
produce a parameter-sensitivity right-hand-side by itself — the field would need
its own sensitivities. Generating $\mathbf{R}^{(p)}$ requires a **parameterized
formulation/material backend** whose residual can be evaluated hypercomplexly (as
here). See [../../docs/output_contract.md](../../docs/output_contract.md).
