# Python residual template

A 1-DOF cubic spring `R = k·u³ − f`, sensitivities of `u` w.r.t. `k` and `f`.

## Edit
- `user_residual.py` — put your residual in `residual(u, params, state, time)`
  (and optionally `tangent(...)`). Use ordinary arithmetic so OTILib can carry
  the derivatives.
- `resasm.yml` — set `problem.unknowns`, your `parameters`, and `sensitivity.order`.
- `solution.npy` — your converged solution vector.

## Run
```
resasm check resasm.yml
resasm run   resasm.yml
```

## Expect
`resasm_output/public/summary.md` with the parameter ranking, and
`resasm_output/private/` with the full `R^(p)` / `U^(p)` arrays. For this demo the
order-1 sensitivities are `du/dk = −1/3` and `du/df = 1/24`.

Requires OTILib (see QUICKSTART_USER.md). Order ≥ 2 needs OTILib; the finite-
difference cross-check confirms order 1.
