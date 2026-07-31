# Example — black-box residual (no linking)

The partner exposes only an **executable** ([mock_solver.py](mock_solver.py))
that speaks the kit's JSON contract. The kit writes `request.json`, runs the
solver, and reads `response.json` — the model stays entirely behind the
executable boundary. The kit performs no overloading here; the executable returns
`R^(1)` itself.

## Run

```bash
python -m partner_kit.python.partner_cli --config config.json
```

## Expected output

```
provider        : blackbox (BlackBoxResidualProvider)
parameters      : ['k']
tangent source  : blackbox-response
R^(1) shape     : (1, 1)
hypercomplex ok : True
  dU/d(k         ) = -3.333333e-01
validation      : True
```

Local validation for a black-box provider:

```
| check | value | tol | result | detail |
| real-residual     | — | 1e-08  | pending | black-box: check inside the executable |
| tangent-vs-FD     | — | 0.0001 | pending | unavailable (black-box or no tangent)  |
| rhs-vs-FD         | — | 0.0001 | pending | black-box: check inside the executable |
| sensitivity-solve | 0 | 1e-08  | PASS    | T dU/dp + dR/dp ~ 0 |
```

The finite-difference checks are **pending** because the kit cannot evaluate the
partner's residual locally — the executable should self-validate (the reference
solver in `templates/blackbox_executable` shows how). The `sensitivity-solve`
check still runs from the returned `T` and `R^(1)`.

## Contract

- Request: `resasm-partner-request/1` (`order`, `solution`, `parameters`,
  `seed_directions`, `time`, `dtime`).
- Response: `resasm-partner-response/1` (`residual_coefficients` = `R^(p)`,
  optional `residual_real`, `tangent`, `diagnostics`).

See [../../docs/residual_provider_contract.md](../../docs/residual_provider_contract.md).
