# Privacy Model

The kit is designed so a collaborator can benefit **even if they share nothing
back except final sensitivity results** — and can share nothing at all while
still using it locally.

## Data that can remain private (never required to leave your machine)

- **source code** (your provider implementation)
- **mesh / geometry**
- **input deck**
- **material model** (UMAT/UEL/constitutive code)
- **state variables** (STATEV, history)
- **the full residual vector** $\mathbf{R}$
- **the tangent matrix** $\mathbf{T}$

The kit only ever *calls* your provider (or your executable). It does not read,
transmit, or persist any of the above outside your machine. The private
`sensitivity_package/` (full arrays) is written locally and is yours to keep.

## Data you may optionally choose to share

From `public_report/` only:

- **validation report** (`validation_summary.md`) — pass/fail of the local checks
- **timing report** (`timing_summary.json`)
- **derivative norms** (`sensitivity_norms.csv`) — $\lVert dU/da_i\rVert$ per parameter
- **selected output sensitivities** (opt-in only, `share_full_arrays: true`)
- **anonymized parameter ranking** (`parameter_ranking.csv`)

The public report is built to contain **no mesh, no source, and no full
residual/tangent**. Parameter names can be anonymized before sharing (rename them
in your config, e.g. `p1, p2, ...`).

## Boundary summary

| Artifact | Where it lives | Shareable? |
|---|---|---|
| provider source / mesh / material | your machine | **no** (never requested) |
| `sensitivity_package/*.npz` (R, T, R^(p), U^(p)) | your machine | your choice (default: keep) |
| `sensitivity_package/metadata.json`, `diagnostics.json` | your machine | your choice |
| `public_report/validation_summary.md` | your machine | yes (no proprietary data) |
| `public_report/parameter_ranking.csv`, `sensitivity_norms.csv` | your machine | yes (norms only) |
| `public_report/timing_summary.json`, `errors.json` | your machine | yes |

## Opt-in sharing of arrays

If — and only if — you set `"share_full_arrays": true` in the config, the public
report additionally includes `shared_sensitivities.npz`. This is off by default.

## What we can and cannot infer

- From norms/rankings we learn **which parameters matter most**, not the model.
- We never receive the mesh, the constitutive law, or the residual/tangent unless
  you deliberately include them.
- The kit makes **no network calls**. It reads your config, runs your provider,
  and writes local files.
