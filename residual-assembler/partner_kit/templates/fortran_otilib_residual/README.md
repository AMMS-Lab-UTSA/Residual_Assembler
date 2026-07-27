# Template — Fortran OTILib residual

Compute arbitrary-order residual sensitivities from a Fortran residual using
OTILib. Choose the integration path that fits your codebase; your source stays
private.

## Options

1. **OTILib Fortran interface** (if your OTILib build provides one)
   Use OTILib's derived type + operator overloads directly: rewrite the residual
   against `type(oti)` and read coefficients along each order-`p` direction.

2. **Source transformation** (Tapenade, etc.)
   Generate higher-order tangent routines. Robust but needs a build step and
   careful handling of the order-by-order state injection.

3. **Generated wrapper / C interop**
   Keep your real Fortran residual; wrap it so a small C/C++ driver built against
   OTILib (see `../cpp_otilib_residual`) evaluates it with OTI scalars via
   `iso_c_binding`. Requires the residual to be callable with a generic scalar
   (or duplicated for the OTI type).

4. **Black-box executable**
   Handle OTI internally and emit `residual_coefficients_by_order` in
   `response.json` (schema in `../../docs/otilib_partner_quickstart.md`). Your
   code never leaves the process.

## Seeding & extraction (same conventions as everywhere)

- Seed all parameters simultaneously: `a_i* = a_i + e_i`.
- Configure `m` bases, `nt = q` truncation order.
- Direction = exponent multi-index `kappa`; true derivative =
  `product(factorial(kappa_i)) * coeff`.

## Order-by-order

For `p > 1`, inject solved `U^(k)` (k<p) into the hypercomplex `u*` before
evaluating order `p`. This is essential — a single final evaluation is not enough.

## History-dependent materials

For plasticity/damage/viscoelasticity/crystal-plasticity you must replay the load
history with the parameters seeded from the first increment (or propagate state
sensitivities). See [../../docs/history_dependent_models.md](../../docs/history_dependent_models.md).

## Reference

Install OTILib from <https://github.com/mauriaristi/otilib.git> (branch `master`,
**GPLv3**, Windows = WSL only); it ships a Fortran (F95+) static-dense
implementation. A minimal first-order dual `type(dual)` (educational, NOT the
production OTI) is in
[../fortran_uel_wrapper/dual_mod.f90](../fortran_uel_wrapper/dual_mod.f90); use it
only to understand the seeding idea — OTILib is the production backend.
