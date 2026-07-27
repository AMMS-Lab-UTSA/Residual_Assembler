# Connecting to the OTIS UMAT source-transformation tool

This framework assembles the residual `R(u, a)` and its parameter derivatives and
solves for sensitivities. Its one material-specific input is `d sigma / d a` — the
stress sensitivity at fixed strain. The **UMAT source transformation** project
("OTIS") produces exactly that field from a real UMAT. This note documents the
interface between the two so they compose without either side reading the other's
source.

See also the methodology overview: [../METHODOLOGY.md](../METHODOLOGY.md) Section 6.

## Why the coupling is needed

`R(u, a)` depends on the parameters only through the stress:

```
R^(p) = d R / d a = integral_Omega  B^T (d sigma / d a)  dV
```

The stress exported from a nominal Abaqus run is frozen at `a0`; it carries no
dependence on `a`. So `d sigma / d a` cannot be recovered from a nominal export — the
material must be evaluated at a perturbed parameter. OTIS does this by transforming
the UMAT so it computes the derivative in hypercomplex (OTI) arithmetic, at machine
precision, in a single evaluation.

## What OTIS produces

Input: a real UMAT plus a compact JSON **transformation contract**. Output: a
`HYPAD_UMAT` (`combined_oti_user.f90` = transformed UMAT + bundled Fortran OTI
library) plus a transformation report and validation results.

The transform, per the contract:

- **promote** selected real variables to `TYPE(ONUMM<m>N<n>)` (an OTI number with `m`
  imaginary directions truncated at order `n`);
- **seed** the chosen independent variable with an OTI direction via `SETIM`;
- **replace** the hand-written Jacobian block with an **extraction** of the OTI
  imaginary coefficients via `GETIM`.

Contract fields (compact form): `case_name`, `source.file`,
`jacobian.independent/dependent/target`, `variables.seed/promote/constant/real`,
`replace`, `otis.ntens/order`.

Two seedings matter to us:

| Seed (`jacobian.independent`) | Extracted derivative | Use in the residual method |
|---|---|---|
| `DSTRAN` | `DDSDDE = d sigma / d epsilon` | the material part of the tangent `K` |
| a material parameter `a` | `d sigma / d a` | the residual sensitivity `R^(p)` |

## The interface: SDV fields in the ODB

The transformed UMAT writes its OTI-extracted derivatives into **state variables
(SDVs)**. The validated `elastic` case writes its 16 `DDSDDE` components into
`SDV2..17` (see `elastic/otis_results.json` in the OTIS project). Abaqus writes SDVs
to the ODB as field output, and this framework's extractor already reads them:

```
scripts/extract_odb_fields.py  ->  fields.json  ->  { "statev": { "<eid>": [ ... per IP ... ] }, ... }
```

So the coupling requires **no new data channel**. The transformed UMAT emits the
derivative through the same SDV -> ODB -> `fields.json` path the framework uses for
stress and state today.

### End-to-end data flow

```
UMAT + contract --(OTIS transform)--> HYPAD UMAT (combined_oti_user.f90)
   -> abaqus job=... user=combined_oti_user.f90   (Abaqus runs the transformed UMAT)
   -> ODB carries S, U, RF, and SDV = d sigma / d a  (and DDSDDE if DSTRAN also seeded)
   -> extract_odb_fields.py  -> fields.json
   -> assemble R^(p) = integral B^T (d sigma / d a) dV        [residual_core, same as the r assembly]
   -> build K (from DDSDDE) ; solve K U^(p) = -R^(p)  -> du/da
   -> chain rule -> d(output)/da
```

One transformed UMAT, seeded with **both** `DSTRAN` (tangent) and the **parameters**
(sensitivity), supplies every material-dependent quantity the method needs.

### The assembler side is built and verified

The consuming side — read `S` and the SDV derivative fields, build `K` and `R^(p)`,
solve, chain-rule — is implemented and verified in
[../examples/residual_sensitivity_c3d8/sensitivity_from_fields.py](../examples/residual_sensitivity_c3d8/sensitivity_from_fields.py).
It evaluates no material model; it only consumes the exported fields, and reproduces
the Abaqus-checked sensitivities to rel `0.0`. The concrete per-IP SDV interface it
expects (a single first-order parameter, C3D8):

```
SDV[0:36]   DDSDDE  (6x6 row-major, Abaqus Voigt 11,22,33,12,13,23)  -> K
SDV[36:42]  d sigma / d a  (Voigt)                                   -> R^(p)
```

A real transformed UMAT must write these slots (the layout is a convention the
contract fixes). Everything upstream of "write to SDV" is the transformer's job;
everything downstream is done and tested.

**What remains on the transformer side:** the current tool seeds `DSTRAN` and writes
`DDSDDE` only (see `src/umat_oti/oti/seed.py`, `extract.py` — both hardwired to
`DSTRAN`/`DDSDDE`). A **parameter seed** (`SETIM` on a material parameter, `GETIM` of
`d sigma / d a` into an SDV) is the addition needed; the OTI machinery it would reuse
already exists. Parameter seeding has been done by hand (the project's poster/deck)
and at the flow-rule level (`OTI_computeflowrule/`), so it is proven, not yet
automated in the packaged tool.

## Field/component conventions to pin when wiring a real case

- **SDV layout**: which SDV slots hold `d sigma / d a` for each parameter, and (if
  seeded) `DDSDDE`. The framework reads them positionally, so the contract's SDV
  mapping must be recorded alongside `fields.json`.
- **Voigt order**: Abaqus `S`/`DDSDDE` order is `(11,22,33,12,13,23)`; the assembler
  uses the same. `d sigma / d a` must be written in that order.
- **Integration-point ordering**: for a non-uniform field the C3D8 IP order
  (`ABAQUS_C3D8_GAUSS`) must match Abaqus' ODB order — verify against a
  spatially-varying single-element job before trusting a non-uniform state. (Uniform
  fields are ordering-blind, which is why the first checks use them.)

## Current state of the connection

- **Verified independently:** the residual assembly (offline machine precision; live
  Abaqus reaction match 7.8e-17), the sensitivity engine (vs Abaqus FD, 1e-5), the
  OTI engine (WSL), and the OTIS transformer (5 UMATs Abaqus-validated on `DDSDDE`).
- **Pending:** (1) a **parameter-seed** transform (the validated OTIS cases seed
  `DSTRAN`; a parameter seed is the same tooling, different contract, not yet
  validated); (2) a **local build** of the OTIS Fortran — Abaqus 2024 here is
  configured for classic `ifort` and only `ifx` (oneAPI 2026) is installed, and the
  non-interactive VS/oneAPI environment setup is unreliable. The OTIS toolchain builds
  and runs on the Linux HPC (`intel/oneapi/2024.2`); the natural path is to run the
  transformed UMAT there, export the ODB, and do the assembly/solve/verify here.
