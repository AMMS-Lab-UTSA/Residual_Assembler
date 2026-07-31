# One-file Abaqus verification — C3D8, standard elastic material

This is the smallest end-to-end check against a **live Abaqus run**: one C3D8
cube, isotropic elastic, pulled uniaxially. Abaqus computes the stress and the
reaction forces; the residual assembler rebuilds the internal nodal forces from
Abaqus' own exported stress and must reproduce the reactions.

You describe the whole check in **one JSON file** and run **one command**. You do
**not** hand-write a mesh file and a separate field file — the tool exports the
fields from the ODB for you.

## The one file you write — [job.json](job.json)

```json
{
  "name": "elastic_c3d8_stress_driven",
  "model": "elastic_c3d8.inp",     // your Abaqus model (mesh + BCs)
  "odb":   "elastic_c3d8.odb",     // the ODB from the solved job
  "mode":  "stress-driven",
  "compare_reactions": true,
  "tol": 1e-6
}
```

| key | meaning |
|---|---|
| `model` | the Abaqus `.inp` (or a neutral `.json`). Gives the mesh and boundary conditions. |
| `odb` | the solved `.odb`. The tool runs `abaqus python` once to export stress `S`, reactions `RF`, and displacements `U`. |
| `fields` | *use instead of `odb`* to point at an already-exported `fields.json` (no Abaqus needed — the offline path). |
| `mode` | assembly mode; `stress-driven` rebuilds the residual from the exported stress. |
| `compare_reactions` | compare the assembled reaction to Abaqus `RF`. |
| `tol` | relative tolerance for the reaction match. |

## How to run it

```bash
# 1. solve the model in Abaqus (any machine with Abaqus) -> elastic_c3d8.odb
abaqus job=elastic_c3d8 input=elastic_c3d8.inp double interactive

# 2. one command does the rest: export fields, assemble, compare to reactions
resasm verify-job job.json
```

That is the whole workflow. `resasm verify-job` reads `job.json`, exports the ODB
fields into `out/fields.json`, assembles the residual, splits free vs. prescribed
DOFs, and compares to the reactions.

## Expected output

```
residual verification: elastic_c3d8_stress_driven  (mode=stress-driven)
  ndof                 = 24
  ||R_free||           = 1.24e-14   (expect ~0 at equilibrium)
  free-residual check  = PASS
  ||reaction (at BC)|| = 1.48e+02
  sign convention      = R = +RF
  rel |R - RF|         = 7.82e-17
  reaction match       = PASS (tol=1e-06)
  overall              = PASS
```

The assembled internal force reproduces Abaqus' reactions to machine precision
(~1e-16), and the free-DOF residual is ~0 (equilibrium). A full report is written
to `out/verification_report.json`.

## Without Abaqus on this machine (offline)

If you already exported `fields.json` on the solver machine, copy it here and use
[job_offline.json](job_offline.json) (it sets `fields` instead of `odb`):

```bash
resasm verify-job job_offline.json
```

## What this example does and does not prove

- **Does**: the residual-assembly and field-exchange path is correct end-to-end
  against real Abaqus output — free-DOF equilibrium and reaction match.
- **Does not**: this is a *standard elastic* material with a *uniform* stress and
  `nlgeom=NO`. It does not exercise the crystal-plasticity UMAT, a non-uniform
  stress (integration-point ordering), or the finite-strain path. Those are the
  next checks.
