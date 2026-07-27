# The Residual Method — methodology, inputs, outputs, and how it connects to the OTIS UMAT transformer

This document is the complete picture of what this framework is, what you put in,
what you get out, how it works, and how it joins to the **UMAT source
transformation** (OTIS) project. It is written to be read start to finish. For the
one-command quick start see [QUICKSTART_USER.md](QUICKSTART_USER.md); for the honest
per-feature state see [STATUS.md](STATUS.md).

---

## 1. The problem, in one paragraph

An Abaqus user cannot "just write their residual function." Abaqus owns the global
residual `R(u, a)` internally and never exposes it. Yet the residual method — the
efficient way to get parameter sensitivities `du/da` from a single solve — needs
exactly that residual and its derivative with respect to the parameters. So the goal
of this framework is: **reconstruct `R(u, a)` from ingredients the user already has,
differentiate it with respect to design parameters, and return the sensitivities of
whatever outputs the user cares about — without the user rewriting or even opening
the element/material source code.**

The residual, element by element, is the finite-element weak form:

```
R(u, a) = F_internal(u, a, q) - F_external(a, t) + F_constraints(u, t)

F_internal,e = integral_over_element  B^T sigma(u, a, q)  dV
```

`B` is the element strain-displacement operator; `sigma` is the stress the material
model produces; `a` are the design parameters; `q` is the internal state (STATEV).

---

## 2. What you put in (inputs)

There are two "levels" of input: the **conceptual ingredients** the method needs,
and the **concrete files** you actually hand over.

### 2.1 The conceptual ingredients

| Ingredient | What it is | Who normally produces it |
|---|---|---|
| **Mesh + boundary conditions** | nodes, connectivity, element types, `*Boundary`, loads | your Abaqus `.inp` |
| **Converged solution `u`** | the displacement field at the converged increment | your Abaqus solve (ODB) |
| **Integration-point stress `sigma`** | Cauchy stress at each IP | your Abaqus solve (ODB); or a UMAT |
| **State / history `q`** | STATEV (history-dependent materials) | ODB |
| **The stress sensitivity `d sigma / da`** | how the stress changes with each parameter, at fixed strain | the **OTIS-transformed UMAT** (or closed form for simple materials) |
| **The tangent `K = dR/du`** | the consistent stiffness | assembled from `DDSDDE`, exported from Abaqus, or built here |
| **Parameters `a`** | the design variables to differentiate against | you name them |
| **Outputs** | the responses you want derivatives of (e.g. `sigma11`, von Mises, a displacement) | you name them |
| **Order** | how many derivative orders (1, 2, ...) | you choose |

A crucial distinction that shapes everything below:

- **`R(u, a)` itself and the tangent `K`** can be *assembled from exported fields* or
  *exported directly from Abaqus*. Neither needs the material re-run.
- **`d sigma / da`** is the one ingredient that **cannot** come from a nominal Abaqus
  run — the exported stress is frozen at the nominal parameter value. To get it you
  must evaluate the material at a perturbed parameter. That is precisely the job of
  the OTIS-transformed UMAT (Section 6).

### 2.2 The concrete files

For **assembling and verifying `R`** (works today, see Section 5):

| File | Role | How you get it |
|---|---|---|
| `model.inp` | mesh, BCs, loads | you already have it |
| `model.odb` | converged `U`, `RF`, `S`, `SDV` | your Abaqus run |
| `fields.json` | the exported fields, solver-neutral | `abaqus python scripts/extract_odb_fields.py --odb model.odb --out fields.json` (done for you) |
| `job.json` | **one file** naming model + odb/fields + mode + tolerance | you write it (see below) |

The single **job file** you write is the whole user-facing input for verification:

```json
{
  "name": "my_check",
  "model": "my_model.inp",
  "odb":   "my_model.odb",
  "mode":  "stress-driven",
  "compare_reactions": true,
  "tol": 1e-6
}
```

Schema: [schemas/resasm_verification_job_v1.schema.json](schemas/resasm_verification_job_v1.schema.json).
Worked example: [examples/abaqus_elastic_c3d8/](examples/abaqus_elastic_c3d8/).

For **parameter sensitivities** you additionally provide the material as an
**OTIS-transformed UMAT** plus a small **transformation contract** naming the
parameters to seed and the order (Section 6).

---

## 3. What you get out (outputs)

### 3.1 The conceptual outputs

| Output | Meaning |
|---|---|
| **`R(u, a)`** | the global residual vector, split into free DOFs (~0 at equilibrium) and prescribed DOFs (equal to the reactions) |
| **`du/da`** | the first-order solution sensitivity to each parameter, from the residual solve `K U^(1) = -R^(1)` |
| **higher orders** | `d^2 u / da^2`, mixed sensitivities — the OTI machinery is arbitrary order |
| **output sensitivities** | `d(sigma11)/da`, `d(vonMises)/da`, `d(displacement)/da`, ... via the chain rule |
| **a verification report** | norms, reaction match, pass/fail — shareable without the model |

### 3.2 The concrete files

| File | What is in it |
|---|---|
| `out/fields.json` | the fields extracted from the ODB (stress, reactions, displacements, SDV) |
| `out/verification_report.json` | `||R_free||`, reaction match, sign convention, pass/fail, tolerances |
| `R.npy` | the assembled residual vector (from `resasm assemble --out R.npy`) |
| sensitivity package (`.json/.npz/.md`) | parameters, derivative orders solved, `du/da`, validation summary (`resasm sensitivity --out`) |

The privacy property: the outputs that leave your machine are **norms and rankings**,
not the mesh or the model. Your model stays local.

---

## 4. How it works (the methodology)

The residual method turns one nonlinear solve into arbitrary-order parameter
sensitivities. The steps:

1. **Solve once, nominally.** Abaqus solves `R(u, a0) = 0` at the nominal parameters
   `a0`. This is the one expensive step, and it is the solve you were already going to
   run. The converged `u`, stress `S`, and state `SDV` are exported.

2. **Assemble the residual `R(u, a0)`.** Using the mesh (`B` operators) and the
   exported stress, the framework rebuilds `F_internal = integral B^T sigma dV` element
   by element and subtracts the external load. On the free DOFs this is ~0
   (equilibrium); on the prescribed DOFs it equals the reactions. This is the
   `stress-driven` assembly path, and it is verified against Abaqus to machine
   precision (Section 5).

3. **Form the parameter derivative `R^(p)`.** The residual depends on the parameters
   only through the stress: `dR/da = integral B^T (d sigma / da) dV`. The stress
   sensitivity `d sigma / da` (at fixed strain) is the material's contribution — for a
   real UMAT it comes from the OTIS-transformed UMAT (Section 6); for elasticity it is
   `sigma / E` in closed form. Higher orders use OTI/HYPAD arithmetic so that a single
   evaluation yields all requested orders at machine precision.

4. **Solve the sensitivity system.** With the tangent `K = dR/du`:

   ```
   K U^(p) = -R^(p)
   ```

   The tangent is formed once and reused; each extra parameter or order is one more
   right-hand side and a back substitution — **no new nonlinear solve, no Abaqus
   re-run.** `U^(1) = du/da`.

5. **Chain-rule to the outputs.** For a user output `g(u, a)`:

   ```
   d g / da = (dg/du) . (du/da) + (dg/da)
   ```

   For a displacement output this is just `du/da`. For a stress output like `sigma11`
   or von Mises, `dg/du` is the material tangent times `B`, and `dg/da` is the direct
   stress sensitivity — both already in hand.

The engine that does steps 3–5 is
[examples/residual_sensitivity_c3d8/sensitivity_engine.py](examples/residual_sensitivity_c3d8/sensitivity_engine.py).
It is material-agnostic: the only material-specific input is `d sigma / da`.

### Why the tangent `K` is not a blocker

`K` is either assembled here from the material tangent `DDSDDE`, or **exported from
Abaqus** (Abaqus can write the element/global stiffness). So of the three ingredients
the method needs, two (`R` and `K`) are assemblable-or-exportable and one
(`d sigma / da`) is the genuinely new piece the OTIS transform provides.

---

## 5. What we have done to test the methodology

Nothing here is claimed because "the code runs." Each claim is a measured number.

### 5.1 Assembly of `R` — verified offline to machine precision

On the C3D8 kernel (`residual_core/formulations/c3d8_kernel.py`), against answers
known in closed form:

| Check | State fed in | Compared against | Rel. error |
|---|---|---|---|
| Divergence-theorem patch | uniform stress, one hex | surface-traction nodal forces | 5e-16 |
| Linear-stress + body load | linear stress field | exact body-force balance | 1e-13 |
| Self-equilibrium | uniform stress | zero net force and moment | 1e-13 |
| Small-strain tangent | perturbation of `u` | finite difference of `F_int` | 1.6e-16 |
| Frame objectivity | rigid rotation | rotated stress, rotated forces | 5e-16 |
| Generic core vs kernel | identical inputs | the original hand kernel, bitwise | 0 |

### 5.2 Assembly of `R` — verified against a live Abaqus run

A single elastic C3D8, 1/8-symmetry, uniaxial. Abaqus computes the stress and
reactions; the assembler rebuilds the internal force from Abaqus' own exported
stress and reproduces the reactions:

```
||R_free||     = 1.24e-14   free-DOF equilibrium   -> PASS
rel |R - RF|   = 7.82e-17   reaction vs Abaqus RF  -> PASS  (R = +RF)
```

Reproduce: [examples/abaqus_elastic_c3d8/](examples/abaqus_elastic_c3d8/) with
`resasm verify-job job.json`.

### 5.3 The sensitivity method — verified against Abaqus finite differences

The residual-method sensitivities matched Abaqus (solved at `E-` and `E+`, central
differenced) for both output kinds:

| Case | Output sensitivity | Engine (exact) | Abaqus FD | rel. err |
|---|---|---|---|---|
| load control | `du/dE` (loaded node) | −2.267574e-09 | −2.267535e-09 | 1.7e-05 |
| load control | `du/dE` (Poisson lateral) | +6.802721e-10 | +6.802674e-10 | 6.9e-06 |
| displacement control | `dσ11/dE` (all 8 IPs) | +1.000000e-03 | +1.000032e-03 | 3.2e-05 |

The ~1e-5 gap is the finite-difference truncation on Abaqus' side; the engine returns
the exact analytic derivative. Reproduce:
[examples/residual_sensitivity_c3d8/](examples/residual_sensitivity_c3d8/).

### 5.4 The OTI/HYPAD engine — verified in WSL

The hypercomplex backend (OTILib) passes its suite in WSL (11 passed with
`RUN_OTILIB_TESTS=1`, so a skip would have failed), including the order-2
recovery-factor proof (`d^2 u / dk^2` recovered from the raw OTI coefficient via the
factorial recovery factor). See [STATUS.md](STATUS.md) Part I-c.

### 5.5 The OTIS UMAT transformer — validated on real UMATs

The companion source-transformation tool has 5 UMATs Abaqus-validated on the
`DDSDDE` path (`elastic`, `code_exp`, `code_imp`, `spin_elas_def`, `spin_elastic`),
reaching 18/19 across the batch. See Section 6.

---

## 6. How this connects to the OTIS UMAT source-transformation tool

The residual method's one missing ingredient is `d sigma / da`. That is exactly what
the **UMAT source transformation** project (the "OTIS" tool) produces.

### 6.1 What the OTIS tool does

It takes a real UMAT plus a small JSON **transformation contract** and emits a
**HYPAD UMAT**: the real constitutive law is preserved, but selected variables are
promoted to OTI numbers (`TYPE(ONUMM<m>N<n>)`), the chosen input is seeded with an
independent OTI direction (`SETIM`), and the derivative coefficients are extracted
(`GETIM`). The math is unchanged; the derivatives fall out at machine precision, any
order, in one evaluation.

```
Real UMAT  +  contract   --(OTIS transform)-->  HYPAD UMAT  +  OTI support library
      seed = DSTRAN   -> DDSDDE = d sigma / d epsilon   (the tangent K's material part)
      seed = a        -> d sigma / d a                  (the residual method's R^(p) ingredient)
```

The contract names the source, the variable to seed (`jacobian.independent`), the
output to read (`jacobian.dependent`), the variables to promote, the lines to
replace, and the OTI order.

### 6.2 The clean seam: SDV fields

The transformed UMAT writes its OTI-extracted derivatives into **state variables
(SDVs)** — the validated `elastic` case writes its 16 `DDSDDE` components into
`SDV2..17`. Abaqus outputs SDVs to the ODB, and **this framework's extractor already
reads SDVs** (`scripts/extract_odb_fields.py` collects `statev`). So:

```
OTIS HYPAD UMAT  writes d sigma / d a into SDVs
      -> Abaqus ODB
      -> extract_odb_fields.py reads the SDV field
      -> assemble  R^(p) = integral B^T (d sigma / d a) dV       [SAME machinery as the r assembly we verified]
      -> solve  K U^(p) = -R^(p)  ->  du/da
      -> chain rule -> output sensitivities
```

Assembling `R^(p)` is the **same** operation as assembling `R` — it just reads the
derivative field instead of the stress field. One OTIS UMAT, seeded with **both**
`DSTRAN` directions (giving `DDSDDE`, hence the tangent `K`) and **parameter**
directions (giving `d sigma / da`, hence `R^(p)`), supplies everything the method
needs, all as exported ODB fields.

### 6.3 What is done and what remains on the connection

- **Done:** the transformer works and 5 UMATs are Abaqus-validated on the `DDSDDE`
  seed; the residual assembler and the sensitivity engine are verified (Section 5);
  the SDV seam is the natural interface both sides already speak.
- **Remaining:** the validated OTIS cases seed `DSTRAN` (giving the tangent). A
  **parameter seed** (giving `d sigma / da`) is the same tooling with a different
  contract, but is not yet among the validated cases. And the local Windows build of
  the OTIS Fortran is blocked (Abaqus is configured for classic `ifort`; only `ifx`
  is present) — the OTIS toolchain builds on the Linux HPC. Details:
  [docs/otis_umat_connection.md](docs/otis_umat_connection.md).

---

## 7. Status in one table

| Capability | State |
|---|---|
| Assemble `R(u, a)` from ingredients (C3D8) | **verified** offline (machine precision) and against a live Abaqus run (7.8e-17) |
| Single-file verification (`resasm verify-job`) | **works** ([examples/abaqus_elastic_c3d8](examples/abaqus_elastic_c3d8/)) |
| Sensitivity method (`R^(p)` -> solve -> chain rule) | **verified vs Abaqus FD** for the elastic case (1e-5) |
| `d sigma / da` for a real UMAT (OTIS transform) | transformer validated on `DDSDDE`; parameter-seed transform + local build pending |
| Tangent `K` (small strain) | assembled + FD-verified; finite-strain mapping pending |
| Higher-order / OTI engine | verified in WSL |

See [STATUS.md](STATUS.md) for the full, honest accounting including limits (element
coverage, load types, MPCs, finite-strain tangent, IP ordering).
