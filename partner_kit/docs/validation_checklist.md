# Validation Checklist

All checks run **locally** on the partner machine (`python/validators.py`). They
need no shared data. Run them via `partner_cli` (they populate
`validation_report.md` / `validation_summary.md`) or call `validators` directly.

## 1. Real residual — `R_free ≈ 0`

At a converged solution the free-DOF residual must vanish.

- Check: `check_real_residual(provider, U, ...)` → `||R_free|| < tol` (default 1e-8).
- Fails if `U` is not converged or the residual sign/scaling is wrong.

## 2. Tangent — `T ≈ dR/dU`

The tangent must match a directional finite difference of the residual.

- Check: `check_tangent_directional(provider, U, T, ...)` →
  `||T·v − (R(U+hv) − R(U−hv))/2h|| / ||·|| < tol` (default 1e-4) over several
  random directions `v`.
- Catches wrong analytic tangents, missing geometric terms, or a stale `T`.

## 3. RHS — `R^(1) ≈ dR/dp`

The generated sensitivity RHS must match a finite difference in the parameter.

- Check: `check_rhs_fd(provider, U, R1, ...)` → per parameter,
  `||R1[:,i] − (R(p+h) − R(p−h))/2h|| / ||·|| < tol` (default 1e-4).
- Confirms the hypercomplex seeding extracted the right derivative.

## 4. Sensitivity — `T dU/dp = −dR/dp`

The solved system residual must vanish.

- Check: `check_sensitivity_solve(T, R1, U1, free_mask)` →
  `||T_ff U1_free + R1_free|| / ||R1_free|| < tol` (default 1e-8).

## 5. (Optional) Output derivatives vs finite differences

If you define a scalar output `g(U, params)`, compare `dg/dp` from the chain rule
(`dg/dU · dU/dp + ∂g/∂p`) against a finite difference — a strong end-to-end check.

- Helper: `check_output_sensitivity_fd(provider, U, U1, output_fn, ...)`.

## Black-box providers

For a black-box executable the kit cannot evaluate the residual itself, so checks
1–3 are reported **pending** — your executable should self-validate (see the
`blackbox_executable` template, which returns diagnostics). Check 4 still runs from
the returned `T` and `R^(p)`.

## Interpreting results

- `PASS` — within tolerance.
- `FAIL` — investigate (unconverged `U`, wrong tangent, wrong parameter mapping,
  non-scalar-generic residual).
- `pending` — not evaluable locally (black-box) or data not provided.

The `validation_summary.md` in the public report contains only the pass/fail
table — no arrays, mesh, or source.
