# Metric definitions and reconciliation

Every published validation number is an **OTI-vs-FD relative discrepancy** — the
relative disagreement between the OTI derivative and a finite-difference estimate
of the *same* derivative. It is deliberately **not** called a "relative error":
there is no analytical derivative to be the truth here, so neither side is "the
right answer." Where the discrepancy is no larger than the FD reference's own
convergence uncertainty (its Richardson gap), OTI and FD **agree to the
resolution of the FD reference** — and one cannot say which is closer to the true
derivative. A smaller discrepancy is not evidence that OTI is more accurate than
FD.

## Why the FCC self-check reports 2.95e-06 and the slide table reports 9.36e-09

These are **different metrics, on different loading paths, over different
quantities**. Both are correct; they answer different questions. Nothing
disagrees.

| | transformer self-check (`umat_transform.py`) | slide table (`build_error_tables.py`, `prog1_path`) |
|---|---|---|
| **quantity** | `DSIGMA_DP` **and** `DSTATEV_DP` | `DSIGMA_DP` only |
| **aggregation** | **max** component (worst single point) | **RMSE** over path and components |
| **loading path** | 5 increments: elastic → plastic → **unload → reload** | 150 increments: **monotonic** uniaxial ramp to 3% |
| **normalization** | per parameter: max \|converged FD column\| over the path | per parameter: max \|FD\| over the path, floored at 1e-3·global scale |
| **worst value (m6_fcc)** | `2.95e-06` (DSTATEV_DP), `1.91e-06` (DSIGMA_DP) | `9.36e-09` (DSIGMA_DP) |

The gap is the product of three deliberate differences:

1. **max vs RMSE.** The self-check reports the single worst point; the slide
   reports a root-mean-square. Over a path, `max ≥ RMSE` always, and the ratio
   grows when the discrepancy is concentrated at a few increments.
2. **5-point harsh path vs 150-point smooth path.** The self-check path
   *unloads and reloads* — the harshest test of a path-dependent derivative,
   with abrupt regime changes — and samples only 5 increments, so the max lands
   on a transition point. The slide path is a smooth 150-step monotonic ramp,
   where the derivative varies gently and the RMSE averages down.
3. **State vs stress derivative.** `2.95e-06` is the **state** derivative
   `DSTATEV_DP` (slip resistances), which the slide table does not include; the
   self-check's `DSIGMA_DP` figure is `1.91e-06`, and even that is a max over the
   harsh path, not an RMSE over the smooth one.

Restricting to the same quantity (DSIGMA_DP): `1.91e-06` (max, 5-step
unload/reload) vs `9.36e-09` (RMSE, 150-step monotonic). Both are far below the
`1e-4` pass tolerance, and both are within — or comparable to — the FD
reference's own Richardson uncertainty at the selected step (see the per-parameter
`fd_reference_uncertainty` in each JSON). So OTI agrees with FD to the FD
reference's resolution in both; the numbers differ only because one is a
worst-case point on a deliberately harsh path and the other is an averaged figure
on a smooth one.

## Every metric, precisely

### Self-check (Program 1 build) — `umat_transform.py`

Full definitions are embedded in each completed contract under
`validation.metric_definitions`. Summary:

| metric | formula | agg. | path | norm. denominator | pass |
|---|---|---|---|---|---|
| `stress_parity_max_rel` | `max\|σ_OTI−σ_orig\| / max\|σ_orig\|` | max | full 5-step | max \|σ_orig\| (floor 1e-30) | ≤1e-8 |
| `statev_parity_max_rel` | `max\|S_OTI−S_orig\| / max\|S_orig\|` | max | full 5-step | max \|S_orig\| (floor 1e-30) | ≤1e-8 |
| `ddsdde_parity_max_rel` | `max\|DDSDDE_OTI−F*_DSTRAN\| / max\|F*_DSTRAN\|` | max | last increment | max \|FD tangent\| (floor 1e-30) | ≤1e-4 |
| `dsigma_max_rel` | `max_c max\|DSIGMA_DP_OTI[:,c]−F*(c)\| / max\|F*(c)\|` | max, worst param | full 5-step | per-param max \|FD col\| (floor 1e-30) | ≤1e-4 **and FD converged** |
| `dstatev_max_rel` | `max_c max\|DSTATEV_DP_OTI[:,c]−F*(c)\| / max\|F*(c)\|` | max, worst param | full 5-step | per-param max \|FD col\| (floor 1e-30) | ≤1e-4 **and FD converged** |

- **Near-zero reference derivatives.** The denominator is floored at `1e-30` (or,
  in the deck tables, at `1e-3 × global scale`). A near-zero FD column inflates
  the ratio; it is reported honestly, never silently set to `0`. Any `NaN`/`Inf`
  is an unconditional failure.
- **Selected FD step & FD-reference uncertainty.** Reported per parameter:
  `selected_h_rel`, and `fd_reference_uncertainty` = the Richardson gap (or the
  successive-change) at that step. `within_fd_uncertainty` is true when the
  OTI-vs-FD discrepancy ≤ ~3× that gap.

### Slide/table (Program 1 & 2 figures) — `residual-assembler/results/`

Definitions embedded in `results/figures/error_tables.json` (`_meta.metric_definitions`)
and `program1_material_validation.json` (`metric_definition`). Summary:

| script / metric | formula | agg. | path | norm. |
|---|---|---|---|---|
| `prog1_stateless` (elastic) | `max\|DSIGMA_DP_OTI[:,p]−F*(p)\| / max(\|F*(p)\|, 1e-3·gs)` | max | single point | per-param max\|FD\|, floor 1e-3·global |
| `prog1_path` (path-dep.) | `RMSE(DSIGMA_DP_OTI−F*) / max(\|F*\|, 1e-3·gs)` | RMSE | 150-step monotonic uniaxial to 3% | per-param max\|FD\| over path, floor 1e-3·global |
| `prog2_residual` | `\|dσ_vM/dp\|_resid − F*\| / max(\|F*\|, 1e-3·gs)` | scalar | point / last of 120-step | per-param \|FD\|, floor 1e-3·global |
| `program1_material_validation` `worst_rel_rmse` | `RMSE(p·dOTI/dp − p·dFD/dp)/max\|p·dFD/dp\|` | RMSE, worst param & block | 60-step uniaxial ramp to 1% | max \|param-weighted FD\| |

In all of the above, `F*` is the **converged** finite-difference reference from
`umat_oti.validation.fd_reference`: a parameter-scaled centered difference whose
step is selected per parameter from the FD sequence's own two-sided plateau and
checked by Richardson extrapolation — never chosen by agreement with OTI.

## Terminology rule applied throughout

- "OTI-vs-FD relative discrepancy" (or "discrepancy") — used everywhere, because
  the reference is FD, not analytical.
- "error" — reserved for a comparison against an analytical derivative (the
  framework has none, so this term is not used for OTI vs FD).
- "OTI agrees with FD within the FD reference's uncertainty" — stated when the
  discrepancy is ≤ ~3× the Richardson gap. In that regime the residual is
  finite-difference truncation, and no claim is made that OTI is *more accurate*
  than FD.
