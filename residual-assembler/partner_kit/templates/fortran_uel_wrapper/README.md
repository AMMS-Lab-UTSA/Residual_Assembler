# Template — Fortran UEL/UMAT residual wrapper

Fortran has no native operator overloading across a legacy codebase, so pick the
integration path that fits your code. All keep your source private; the kit only
consumes the residual coefficients you produce.

## Options (choose one)

1. **Source transformation** (Tapenade, etc.)
   Generate a differentiated routine `eval_residual_d` that returns `dR/dp` for a
   seeded parameter. Highest fidelity; needs a build step.

2. **Operator-overloaded derived type**
   Define a `type(dual)` with overloaded `+ - * / **` and rewrite the residual
   against it (see `dual_mod.f90` skeleton). The imaginary component carries the
   derivative — the same idea as the C++ `Dual1`.

3. **Generated wrapper**
   Wrap your existing real `UEL`/`UMAT` in a small driver that computes `dR/dp`
   for the requested parameters, then emits `response.json` (schema in
   [../../docs/residual_provider_contract.md](../../docs/residual_provider_contract.md)).

4. **Black-box finite-difference fallback — validation only**
   If you cannot differentiate the routine, compute `dR/dp` by finite differences
   inside your executable. Use this only to *validate* the pathway; it is subject
   to step-size/truncation error and must not be presented as exact.

## What the kit needs from you

Either implement a Python `ElementResidualProvider` that calls your compiled
routine (via f2py/ctypes), or ship a black-box executable that writes
`response.json` with `residual_coefficients` (= `R^(p)`), optional `residual_real`
and optional `tangent`.

## History-dependent materials

For plasticity/damage/viscoelasticity/crystal-plasticity you **must** replay the
load history with the parameter seeded from the first increment (or propagate
state sensitivities). See
[../../docs/history_dependent_models.md](../../docs/history_dependent_models.md).
Final-step-only differentiation is invalid.
