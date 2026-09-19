# Python residual template

This template is the starting point for a direct residual (Path C): you write
`R(u, params)` in Python and Residual_Assembler computes its parameter
sensitivities to any order with OTILib. Copy it with
`resasm init --template python --out my_case`.

The shipped model is a 1-DOF cubic spring `R = k·u³ − f`, with sensitivities of
`u` with respect to `k` and `f`.

## Edit
- `user_residual.py` — put your residual in `residual(u, params, state, time)`
  (and optionally `tangent(...)`). Use ordinary arithmetic so OTILib can carry
  the derivatives.
- `resasm.yml` — set your `parameters` and `sensitivity.order`
  (`problem.unknowns` is optional; it is inferred from the solution vector).
- `solution.npy` — your converged solution vector.

## Run
```
resasm check resasm.yml
resasm run   resasm.yml
```

## Expect
`resasm_output/public/summary.md` with the parameter ranking, and
`resasm_output/private/` with the full `R^(p)` / `U^(p)` arrays. For this demo the
order-1 sensitivities are `du/dk = −1/3` and `du/df = 1/24`, and at order 2
`d²u/dk² = 2/9` (read `U_derivatives`, not the raw coefficients).

Requires a genuine OTILib build (`scripts/setup_otilib.sh` or
`docs/OTILIB_VENV.md` in the Residual_Assembler repository). The RHS
finite-difference cross-check confirms order 1.
