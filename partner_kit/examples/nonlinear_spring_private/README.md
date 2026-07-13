# Example — private spring (global provider)

A partner's proprietary model stand-in: a 1-DOF cubic spring `R = k u^3 - f`,
implemented in [provider.py](provider.py) as a `GlobalResidualProvider`. The
source, values, and residual never leave the machine; the kit only calls the
provider and reads `dR/dk` from the imaginary part of a dual evaluation.

## Run

```bash
python -m partner_kit.python.partner_cli --config config.json
```

## Expected output

```
provider        : private_spring (PrivateSpring)
parameters      : ['k']
tangent source  : provider-get_tangent
R^(1) shape     : (1, 1)
hypercomplex ok : True
  dU/d(k         ) = -3.333333e-01
validation      : True
outputs         : out/{sensitivity_package,public_report}
```

`du/dk = -u/(3k) = -2/(3·2) = -1/3` — matches the analytic value; the local
validation (real residual, tangent-vs-FD, RHS-vs-FD, sensitivity-solve) all pass.

## Outputs

- `out/sensitivity_package/` — **private**: `metadata.json`, `dof_map.json`,
  `parameter_map.json`, `residual_real.npz`, `tangent.npz`, `rhs_order_1.npz`,
  `solution_sensitivities_order_1.npz`, `validation_report.md`, `diagnostics.json`.
- `out/public_report/` — **shareable**: `validation_summary.md`,
  `parameter_ranking.csv`, `sensitivity_norms.csv`, `timing_summary.json`,
  `errors.json` (no mesh, source, or full residual/tangent).
