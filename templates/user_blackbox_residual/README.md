# Black-box executable template

For users who cannot share or link code. Your model stays inside `my_solver.py`
(or any executable). The framework only exchanges files.

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

Request/response contract: see `docs/user_input_contract.md` and
`docs/user_output_contract.md`.
