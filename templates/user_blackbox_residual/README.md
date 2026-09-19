# Black-box executable template

This template is the starting point for a black-box residual (Path B) at
order 1, for users who cannot share or link code. Your model stays inside
`my_solver.py` (or any executable); the framework only exchanges files. Copy it
with `resasm init --template blackbox --out my_case`; for order 2 and above use
the `blackbox-order2` template.

## Edit
- `my_solver.py` — replace `_residual` with your model. Read `request.json`,
  write `response.npz` with `R_order_<p>` arrays (and optionally `tangent`).
- `resasm.yml` — set `problem.unknowns`, `parameters`, `sensitivity.order`, and
  the `residual.command` (use `{request}` / `{response}` placeholders).
- `solution.npy` — your converged solution vector.

## Run
```
resasm check resasm.yml
resasm run   resasm.yml
```

## Expect
`resasm_output/public/summary.md` with the parameter ranking. This path does
**not** require OTILib on the framework side — the executable produces the
`R^(p)` coefficients itself (here via its own finite differences). For order ≥ 2
your executable must also consume `u_star_coefficients` from the request to
rebuild `u*`.

Request/response contract: `docs/user_input_contract.md` and
`docs/residual_provider_contract.md` in the Residual_Assembler repository.
