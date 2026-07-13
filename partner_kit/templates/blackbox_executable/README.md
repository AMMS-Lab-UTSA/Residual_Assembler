# Template — black-box executable

Ship an executable that speaks the kit's JSON contract; keep everything else
private. The kit calls:

```
your_solver --input request.json --output response.json
```

- Read `request.json` (`resasm-partner-request/1`): `order`, `solution`,
  `parameters`, `seed_directions`, `time`, `dtime`.
- Compute `R^(1)` however you like (your own AD, a dual/OTI type, source
  transformation, or a finite-difference fallback *for validation only*).
- Write `response.json` (`resasm-partner-response/1`) with
  `residual_coefficients` (= `R^(p)`, `ndof x m`), optionally `residual_real`,
  `tangent`, and `diagnostics`.

Schemas: [../../docs/residual_provider_contract.md](../../docs/residual_provider_contract.md).

## `reference_solver.py`

A dependency-free reference implementation with a **finite-difference fallback**
(so it works for any residual you drop in) and an optional analytic hook. Adapt
the `residual(u, params)` function to call your model; the request/response
plumbing stays as-is. Register it with:

```json
{ "blackbox": { "command": ["python", "reference_solver.py"],
                "ndof": 1, "parameters": ["k"],
                "parameter_values": {"k": 2.0, "f": 16.0} },
  "parameters": ["k"], "order": 1, "solution": [2.0], "output_dir": "out" }
```

C++ equivalent: [../cpp_global_residual](../cpp_global_residual).
