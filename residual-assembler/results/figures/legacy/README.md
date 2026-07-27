# Legacy validation artifacts (archived)

These are the **original fixed-step** finite-difference results, kept for
provenance. They are **superseded** by the canonical adaptive-FD artifacts in the
parent directory.

| artifact | method | status |
|---|---|---|
| `error_tables_fixed_step_v0.json` | single hard-coded FD step (P1 path 1e-5, P1 point 1e-6, P2 1e-5/1e-6) | **legacy** |
| `../error_tables.json` | `umat_oti.validation.fd_reference` — converged, per-parameter step, Richardson-checked | **canonical** |

Why they differ: a single fixed FD step is not converged for every parameter of
every model (a rate exponent of 20 in the FCC crystal needs a far finer step than
an elastic constant). The canonical method selects the step per parameter from
the FD sequence's own convergence, so the OTI-vs-FD discrepancy it reports is
against a *converged* reference. The headline worst-case values move from
`4.4e-08 / 4.5e-08` (legacy) to smaller, converged-reference figures; both are
OTI-vs-FD relative discrepancies, not errors against an analytical derivative.

Every canonical JSON carries a `_meta`/`provenance` block with `method_version`,
git commits, environment, the FD ladder, tolerances, and metric definitions, so
the two methodologies are never confused.
