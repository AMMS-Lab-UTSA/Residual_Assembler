# Partner Output Contract

The kit writes two trees under `output_dir`: a **private** package (full arrays,
keep local) and a **public** report (shareable, no proprietary data).

## Private — `sensitivity_package/`

```
sensitivity_package/
  metadata.json                         provider, ndof, order, algebra, parameters,
                                        tangent source, system_ready, timings
  dof_map.json                          ndof + free/prescribed mask
  parameter_map.json                    parameter name -> imaginary basis index
  residual_real.npz                     R (order 0), if evaluable locally
  tangent.npz                           dense T ...
  tangent_operator_info.json            ... OR this, when T is matrix-free/absent
  rhs_order_p.npz                       R  = R^(p)   and   rhs = -R^(p)
  solution_sensitivities_order_p.npz    U  = U^(p)   (present only if T available)
  validation_report.md                  full local validation table
  diagnostics.json                      RHS diagnostics, errors, timings
```

This tree may contain enough to reconstruct sensitivities; it is meant to stay on
the partner machine. Sharing it is entirely optional.

## Public — `public_report/`

```
public_report/
  validation_summary.md      pass/fail of local checks (no arrays)
  timing_summary.json        wall-clock per stage
  parameter_ranking.csv      rank, parameter, norm, metric   (most influential first)
  sensitivity_norms.csv      parameter, norm, metric          (||dU/dp|| or RHS proxy)
  errors.json                any non-fatal issues / notes
```

The public report is built to contain **no mesh, source, or full
residual/tangent**. Parameter names may be anonymized (rename in the config).
Array sharing is **opt-in** only: setting `"share_full_arrays": true` adds
`shared_sensitivities.npz`.

## Metric semantics

- If a tangent was available, per-parameter norm is $\lVert dU/da_i\rVert$
  (`metric = ||dU/dp||`).
- If no tangent was available, the kit still exports $\mathbf{R}^{(1)}$ and ranks
  by $\lVert \mathbf{R}^{(1)}_{:,i}\rVert$ (`metric = ||R^(1)_col|| (RHS proxy)`),
  clearly labelled — the true sensitivity solve is deferred until a tangent is
  supplied.

## Relationship to our framework's output contract

This mirrors the core framework's
[`residual_core/docs/output_contract.md`](../../residual_core/docs/output_contract.md)
(same target system $\mathbf{T}\,\mathbf{U}^{(p)} = -\mathbf{R}^{(p)}$, same
serialization idea) but adds the **private/public split** required when the model
and residual must never leave the partner machine.

## Who generates $\mathbf{R}^{(p)}$? (the honest boundary)

- A **stress-driven** style residual (assembled from an exported field) yields
  $\mathbf{R}$ but **cannot** produce a parameter-sensitivity RHS by itself.
- $\mathbf{R}^{(p)}$ requires a **parameterized residual evaluable hypercomplexly**
  — i.e. one of the provider levels here. The kit does not claim to compute
  sensitivities without a locally evaluable residual.
