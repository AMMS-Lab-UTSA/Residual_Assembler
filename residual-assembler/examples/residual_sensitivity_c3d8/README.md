# Residual-method parameter sensitivity on a C3D8 — verified against Abaqus

This is the smallest complete demonstration of the **sensitivity** half of the
framework: from one converged elastic solve it computes parameter sensitivities of
displacement and stress outputs, and confirms them against a live Abaqus run.

It is deliberately built on a material we control (isotropic elastic, where
`d sigma / dE = sigma / E` is exact) so every number is checkable. That analytic
`d sigma / dE` is a **stand-in** for the field an OTIS-transformed UMAT supplies for
a real material — see [../../docs/otis_umat_connection.md](../../docs/otis_umat_connection.md).
The engine does not change when the material changes; only the source of
`d sigma / da` changes.

## What it does (the residual method)

```
assemble  R^(1) = d r / dE = integral B^T (d sigma / dE) dV      (parameter derivative of the residual)
build     K = d r / d u                                          (the tangent)
solve     K U^(1) = -R^(1)                                       ->  du/dE
chain     d g / dE = (dg/du) . (du/dE) + (dg/da)                 ->  d(sigma11)/dE, d(vonMises)/dE
```

One linearization, no re-solve. See
[sensitivity_engine.py](sensitivity_engine.py).

## Run it

Self-check (no Abaqus — engine vs its own finite difference and the analytic answer):

```bash
python -m examples.residual_sensitivity_c3d8.sensitivity_engine
```

Verify against Abaqus (needs Abaqus):

```bash
python -m examples.residual_sensitivity_c3d8.verify_against_abaqus gen      # Case A inps + engine prediction
# run Abaqus at E-,E+ and extract fields with scripts/extract_odb_fields.py,
# saving fieldsA_m.json / fieldsA_p.json into abaqus_runs/
python -m examples.residual_sensitivity_c3d8.verify_against_abaqus compare
# Case B (stress output): genB / compareB
```

## Result (measured)

Two loading cases cover both output kinds. The engine returns the exact analytic
derivative; Abaqus is central-differenced across two solves, so the ~1e-5 gap is
the finite-difference truncation, not engine error.

| Case | Output sensitivity | Engine | Abaqus FD | rel. err |
|------|--------------------|--------|-----------|----------|
| A — load control | `du/dE` (loaded node, x) | −2.267574e-09 | −2.267535e-09 | 1.7e-05 |
| A — load control | `du/dE` (lateral, Poisson) | +6.802721e-10 | +6.802674e-10 | 6.9e-06 |
| B — displacement control | `dσ11/dE` (all 8 IPs) | +1.000000e-03 | +1.000032e-03 | 3.2e-05 |

Case A exercises the sensitivity solve itself (nonzero `du/dE`); Case B exercises
the stress-output chain rule (nonzero `dσ11/dE`). The engine's own
finite-difference self-checks are tighter (1e-8 and 1e-13).

## Consuming a material formulation's outputs (the assembler, not the material)

[sensitivity_from_fields.py](sensitivity_from_fields.py) is the same method, but it
reads the derivatives from **exported ODB/SDV fields** instead of computing
`d sigma / dE` in closed form — i.e. it consumes what an OTI-overloaded UMAT writes,
and evaluates no material model itself. This is the residual assembler proper:

```
overloaded UMAT (seeded strain + parameter) writes per integration point:
    S            = sigma
    SDV[0:36]    = DDSDDE  = d sigma / d epsilon     -> assembler builds K = integral B^T D B dV
    SDV[36:42]   = d sigma / d a                     -> assembler builds R^(p) = integral B^T (d sigma/da) dV
assembler: solve K U^(1) = -R^(1) -> du/da ; chain-rule -> output sensitivities
```

Run it:

```bash
python -m examples.residual_sensitivity_c3d8.sensitivity_from_fields
```

It reproduces the direct engine's `du/dE` and `dσ11/dE` to rel `0.0` — the assembler
forms the residual sensitivities **purely from the exported fields**. The SDV layout
above is the interface a material formulation (an overloaded UMAT) must fill; see
[../../docs/otis_umat_connection.md](../../docs/otis_umat_connection.md).

## Scope

Elastic material, single C3D8, `nlgeom=NO`. It proves the residual-method machinery
end to end against Abaqus. The assembler consumes the material formulation only
through the exported `S` + SDV fields, so any material (any overloaded UMAT that
writes `d sigma / da` to SDVs) plugs in without changing the assembler.
