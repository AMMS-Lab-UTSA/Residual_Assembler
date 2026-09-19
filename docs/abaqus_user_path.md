# The Abaqus user path

This page is for Abaqus users. It explains what Residual_Assembler can do with
an Abaqus model: compute sensitivities of a finished analysis, assemble the
residual from the fields Abaqus exports, and which model features cannot be
reconstructed. For a first run see [QUICKSTART_USER.md](../QUICKSTART_USER.md).

## The premise

**Abaqus does not expose the global residual.** There is no output request, no
ODB field, no user-subroutine hook that hands you `R = F_int - F_ext`. Abaqus
forms it internally, uses it, and discards it. A UMAT does not help: it returns
stress, `DDSDDE` and `STATEV` at a *single material point*, and knows nothing
about shape functions, the B-operator, the integration rule, or how those point
quantities become element nodal forces (`residual_core/docs/limitations.md` §1).

So we do not ask Abaqus for `R`. **We reconstruct `R` from the ingredients Abaqus
*does* export.**

```
    R(u, a)  =  F_internal(u, a, q)  -  F_external(a, t)  +  F_constraints(u, t)
    R_e      =  ∫_Ωe  Bᵀ σ dΩ  -  f_e^ext
```

Every term on the right is recoverable from the `.inp` plus an ODB export — for
the subset of Abaqus models listed under [Limits](#what-cannot-be-reconstructed).
Outside that subset, the reconstruction is **incomplete**, and this tool should
not be trusted to reproduce Abaqus' reactions.

---

## Sensitivities of a finished analysis

If what you want is the parameter sensitivities of an analysis you have already
run, you do not need the rest of this page. Build the material once with the
companion UMAT-OTI (`umat-oti-provider build`), then:

```bash
resasm request --model Analysis.inp --odb Analysis.odb \
    --material OTI_UMAT.obj --request sensitivity_request.json --out results
```

The command exports every ODB increment with Abaqus Python, replays the
compiled material at every integration point and increment, checks the replayed
stress, state and reactions against the ODB, assembles `K` and `dR/dp`, and
solves for the requested sensitivities. It handles small-strain C3D8 analyses
with one static step, one user material, ramped boundaries (including nonzero
prescribed displacements) and concentrated loads; models outside its bounded
single-material engine go to `resasm history`, which says so. See
[REQUEST_INTERFACE.md](REQUEST_INTERFACE.md),
[REPLAY_HISTORY.md](REPLAY_HISTORY.md) and the full-size cantilevers in
[examples/cantilevers](../examples/cantilevers/README.md).

The rest of this page is about assembling the residual itself.

## Minimum you must supply

| ingredient | where it comes from | required for |
|---|---|---|
| **`.inp`** — mesh, connectivity, sections, materials, BCs, `*Cload` | your job file | everything |
| **ODB export of `U`** (nodal solution) | `scripts/extract_odb_fields.py` / `residual_core/io/abaqus_odb_export.py` | current-configuration (nlgeom) assembly; reaction comparison |
| **the material** — either the UMAT source **or** an exported stress field | your subroutine, or the ODB `S` field | `F_internal` |
| **loads / BCs / stimuli** | the `.inp` (parsed) | `F_external`, the free/prescribed partition |
| **state variables** (`SDV` / `STATEV`) | the ODB | history-dependent materials only |

Two ways to get σ, and they are not equivalent:

| σ from | mode | needs a UMAT? | gives a tangent? | gives sensitivities? |
|---|---|---|---|---|
| Abaqus' exported IP stress | `stress-driven` | no | **no** | **no** |
| a material update from `u` | `material-replay` | yes (compiled) | yes | not through the recipe (see [below](#sensitivities-which-route)) |

---

## The concrete commands

### 1. Inspect the model — what can this tool do with your `.inp`?

```bash
resasm inspect      model.inp
resasm inspect      model.inp --detail        # per-element backend selection
resasm requirements model.inp --mode stress-driven
resasm requirements model.inp --mode material-replay
resasm doctor       model.inp                 # per-mode readiness
resasm backends                               # every backend's declared limits
```

Real output for the shipped crystal-plasticity example
(`sources/permissive/ngrilli_Oxford_Crystal_Plasticity/ExampleInputFiles/HCPnoTwin/Compression111.inp`,
216 nodes, 125 C3D8, `*User Material` CPuranium, 11 constants, `*Depvar 125`):

```
Elements:
  - C3D8 elements (x125): supported by solid_c3d8_finite_strain backend
Materials:
  - CPuranium: UMAT adapter required; source file not provided
Required user inputs:
  1. Provide UMAT source or an exported field for material 'CPuranium'
     (or run stress-driven mode).
```

`resasm requirements --mode <mode>` reports the **single minimum missing item**,
not a checklist.

### 2. Run the job and export the fields (needs Abaqus)

```bash
python scripts/run_abaqus_validation.py --job my_job --inp model.inp --user umat.for
abaqus python scripts/extract_odb_fields.py --odb my_job.odb --out fields.json
```

`scripts/extract_odb_fields.py` writes the **assembler field format**:

```json
{ "stress_ip":     { "<eid>":  [[s11,s22,s33,s12,s13,s23], "… one row per IP"] },
  "reactions":     { "<node>": [rf1, rf2, rf3] },
  "displacements": { "<node>": [u1, u2, u3] },
  "statev":        { "<eid>":  [["… sdv"], "… per IP"] },
  "meta":          { "odb": "…", "step": "…", "frame": 0 } }
```

There is a **second, different** exporter:
`residual_core/io/abaqus_odb_export.py` (run as `abaqus python
residual_core/io/abaqus_odb_export.py -- --odb job.odb --out fields.json
--frames all`). It writes the CONTRACT §5 format — `{nodes, elements, frames:[{U,
RF, S, SDV}]}`, multi-frame — which is what `residual_core/stress_driven_residual.py`
consumes. **The two formats are not interchangeable.** Match the exporter to the
consumer:

| exporter | format | consumer |
|---|---|---|
| `scripts/extract_odb_fields.py` | flat `stress_ip` / `reactions` / `displacements` / `statev` | `resasm assemble --fields`, `ResidualProblem.attach_results`, `scripts/compare_residuals.py` |
| `residual_core/io/abaqus_odb_export.py` | CONTRACT §5, per-frame `U`/`RF`/`S`/`SDV` | `residual_core/stress_driven_residual.py --fields` |

Both are Abaqus-Python (`odbAccess`) scripts. `residual_core/io/abaqus_odb_export.py`
is the exporter `resasm request` runs on real ODBs (tested with Abaqus
2021.HF5); the history engine uses its own `.npz` exporter,
`residual_core/replay/odb_export_npz.py`.

### 3. Assemble the residual offline

```bash
resasm assemble model.inp --mode stress-driven --fields fields.json
resasm assemble model.inp --mode stress-driven --fields fields.json --out R.npy
resasm verify   model.inp --fields fields.json      # ||R_free|| and ||reaction||
```

Or, for the residual-vs-Abaqus comparison with the sign check:

```bash
python scripts/compare_residuals.py --model model.inp --fields fields.json --mode stress-driven
```

`compare_residuals.py` prints `||R_free||` (expect ≈ 0 at a converged frame), the
reaction norm at prescribed DOFs, and — if `reactions` is in the export — the
detected sign convention and `rel |R − RF|` / `rel |R + RF|`.

Or, for the **current-configuration** (nlgeom-consistent) assembly:

```bash
python residual_core/stress_driven_residual.py --inp model.inp \
       --fields fields.json --frame last --mode finite
```

**This distinction matters.** `resasm assemble --mode stress-driven` hard-sets
`options['config'] = 'small'` (`ui/wizard.py::_run`) — the internal force is
integrated in the **reference** configuration and the exported `U` does not enter
it. Abaqus' CP examples are `nlgeom=YES`, and the UMAT returns **Cauchy** stress,
which must be integrated over the **current** configuration. For a finite-strain
job, use `stress_driven_residual.py --mode finite`, which does use the exported
`U`. (At `U = 0` the two paths agree exactly; the difference grows with the
deformation.)

Also note: when `--fields` is a JSON **path**, `attach_results` keeps only
`stress_ip`. Exported `reactions`, `displacements` and `statev` in that same file
are **discarded** by `resasm assemble`. `compare_residuals.py` re-reads
`reactions` from the file itself, which is why it can do the RF check.

### 4. Material replay (needs the compiled UMAT)

```bash
resasm assemble model.inp --mode material-replay --subroutine umat.for --tangent
python scripts/compare_umat_replay.py ...
```

For a `*User Material`, `ResidualProblem` refuses to fabricate a stress:

> `material-replay for user material 'CPuranium' needs the compiled UMAT (Intel
> ifort + Abaqus). Use mode='stress-driven' with an exported field for an offline
> residual.`

The Grilli CP UMAT **does not compile under gfortran** (Cray-pointer/`target` twin
arrays + an ifort `trace()` kind mismatch). It needs **Intel ifort + Abaqus (MKL)**.
With gfortran only a *mock* elastic UMAT runs through this path, so it has
proven its *plumbing* (5-increment `STATEV` history, `/UMPS/` common block
persisted across increments), not crystal plasticity.

A UMAT that UMAT-OTI has transformed and compiled takes the other route: the
analysis replay links the compiled provider and replays it against the ODB at
every integration point and increment (see
[Sensitivities of a finished analysis](#sensitivities-of-a-finished-analysis)).

---

## What CANNOT be reconstructed

If your model contains any of these, the reconstructed `R` is **incomplete** and
will not match Abaqus' reactions. The tool does not warn you — it just assembles
a residual that is missing a term.

| feature | status | consequence |
|---|---|---|
| **`*Dsload`** (pressure / distributed surface load) | parsed into `AbaqusModel.dsloads`, **never assembled** (`core/loads.py` handles `*Cload` only) | `F_external` is missing the surface traction |
| **body forces / gravity** | not parsed, not assembled | `F_external` incomplete |
| **`*Amplitude`** | not a parsed keyword; the BC/load records carry the amplitude *name*, but `Assembler.assemble` calls `external_force(...)` with **no amplitude map** | load and BC scaling are never applied |
| **contact** | not parsed, not assembled | contact forces missing entirely |
| **`*Equation` / MPC / tie** | parsed onto `model.equations`, **never applied** | the constrained DOFs are treated as free |
| **non-C3D8 solids** (`C3D8R`, `C3D20R`, `C3D4`, `C3D10`) | no backend — "planned" | those elements are **silently skipped** (`diag['skipped_no_formulation']`); `resasm inspect` reports "backend missing" |
| **shells** (`S3`/`S4`/`S4R`/…) | `shell_placeholder` — contract only; `eval_element` raises | recognised, not assemblable |
| **user elements (UEL)** | `uel_direct` adapter exists but is registered *unconfigured*; you must supply the Python callable | not assembled until wired |
| **any state Abaqus does not export** | — | a history-dependent replay cannot be reproduced |
| **non-identity `*Instance` transform** | parser raises `NotImplementedError` | refuses to parse rather than silently mis-place the mesh |
| **integration-point ordering** | `ABAQUS_C3D8_GAUSS` uses lexicographic (ξ₁-fastest) order. The internal `sigma_ip[k] ↔ points[k]` pairing is verified offline; the match to **Abaqus' export index** is checked by the replay engines, which compare replayed stress with the ODB at every integration point of non-uniform fields | the stress-driven scripts have no ODB check of their own; compare a spatially varying stress before trusting a non-uniform state through them |

Supported today, end to end: a **single-instance, no-transform, C3D8-only model
with `*Boundary` (Dirichlet/symmetry) and `*Cload`**. That is the honest envelope.
The analysis replay has its own, similar envelope, and refuses anything outside
it by name ([REPLAY_HISTORY.md](REPLAY_HISTORY.md)).

---

## Stress-driven mode is the verification path, not a sensitivity path

Stress-driven mode takes Abaqus' **own** exported integration-point stress and
Abaqus' **own** exported `U`, assembles `∫Bᵀσ dΩ` outside Abaqus, and checks:

- at free DOFs: `R = F_int − F_ext ≈ 0` (Abaqus converged, so this must hold);
- at prescribed DOFs: the assembled internal force `= ± RF` (Abaqus' reaction).

If both hold, **our finite-element assembly reproduces what Abaqus does
internally** — independently of any material model. That is the gold-standard
check, and it is why this mode exists.

It is **not** a sensitivity path. `stress_driven_c3d8` takes σ as a *frozen
exported field*: it does not depend on the material parameters `a`, so
`∂R/∂a = 0` by construction. It also returns no tangent (`k_e = None`;
`TangentSource.UNAVAILABLE`).

### What is verified, and what is not

| claim | state |
|---|---|
| the assembly math (`∫Bᵀσ dΩ`) is correct | ✅ verified **offline to machine precision**: divergence-theorem patch test 5e-16; linear-stress patch + body-load test 1e-13; finite-strain force frame-objective to 5e-16 and reducing exactly to small-strain at `u = 0` |
| the small-strain element tangent `∫BᵀDB` | ✅ FD-verified, relative error 1.6e-16 |
| the `.inp` parse of a real CP job | ✅ verified (216 nodes / 125 C3D8 / 11 constants / 125 SDV, incl. the part/assembly `Set-1` name collision) |
| **the assembled `R` vs Abaqus' reaction forces** | ✅ for the replay engines, which assemble reactions from the replayed stress and compare them with the ODB's `RF` at every increment. `compare_residuals.py` and `stress_driven_residual.py --fields` make the same comparison from exported stress; no run of them against a real ODB is recorded here. See `STATUS.md`. |
| the finite-strain `DDSDDE → AMATRX` mapping | ⏳ open; the finite-strain element tangent is approximate pending the Abaqus comparison |
| the real CP UMAT | ⏳ needs Intel ifort + Abaqus; does not build with gfortran |

Do not read "the math is verified" as "it matches Abaqus". Those are different
statements. The second has been measured for the replay engines, not for the
stress-driven scripts.

---

## Sensitivities: which route

**The assembly recipe (Path A) does not OTI-differentiate a real C3D8 + UMAT
model.**

Why: the recipe's sensitivity engine (`core/oti_rhs_provider.py`) produces
`R^(p)` by evaluating the element residual with hypercomplex (OTI) numbers and
reading the coefficients off. That requires the element and material kernels to
be written in generic arithmetic. The C3D8 kernels are NumPy **float** kernels
(`np.zeros((6,6), dtype=float)` in `core/voigt.py::isotropic_D`, `np.zeros((8,6))`
for the IP stress slab, `np.asarray(dofs, dtype=float)` at every element entry
point, `np.asarray(r_e, float)` in the assembler's scatter). A float array cannot
carry an OTI number: it either raises or truncates to the real part, destroying
the derivative. `Formulation.oti_differentiable` records this: it is `True` for
exactly two proof backends (`nonlinear_spring1`, `nonlinear_bar1`) and `False`
for every solid.

The routes that do give sensitivities:

1. **The analysis replay (`resasm request`, `resasm history`).** UMAT-OTI
   source-transforms the UMAT to OTI arithmetic and compiles it; the replay
   engines link the compiled provider and carry the derivative arrays alongside
   the real ones, so the float kernels never have to hold an OTI number:

   > **OTI-transformed UMAT + this assembler = an OTI-differentiable assembled
   > residual.**

   This is the shipped route for Abaqus analyses; see
   [Sensitivities of a finished analysis](#sensitivities-of-a-finished-analysis).

2. **Bounded finite-strain neo-Hookean** (`resasm sensitivity` with the example
   configuration in [examples/finite_strain_c3d8](../examples/finite_strain_c3d8/README.md)):
   first-order material-parameter sensitivities through a dedicated OTILib
   right-hand side.

3. **Path B, black-box.** Your solver (or your own driver around it) returns the
   order-`p` residual **coefficients**; the framework only solves
   `T U^(p) = −R^(p)`. Your model never leaves your machine. Contract:
   [blackbox_order2_contract.md](blackbox_order2_contract.md), and note the rule:
   return *Taylor coefficients*, not derivatives (at order ≥ 2 they differ by
   `κ!`).

Making the solid kernels themselves OTI-safe would open the recipe path to
C3D8 sensitivities too. It is an engineering task, not a physics limitation,
and has not been done.

---

## Remaining Abaqus comparisons for the assembly framework

Abaqus 2021.HF5 has been used for the replay engines. For the assembly
framework's own scripts, `STATUS.md` and
`residual_core/docs/abaqus_validation_roadmap.md` list what remains:

1. Run a job and export `fields.json` (`--frames all`).
2. Stress-driven check at converged frames with `compare_residuals.py` or
   `stress_driven_residual.py --fields`: expect `||R_free|| ≈ 0` and internal
   force at prescribed DOFs `= ±RF`.
3. Build the Oxford CP UMAT (ifort + Abaqus), replay it increment by increment
   through `material-replay`, compare `STRESS ↔ S` and `STATEV ↔ SDV`, then
   close the loop by feeding the replayed stress back through the (already
   verified) assembler.

---

## See also

- [residual_assembly_recipe.md](residual_assembly_recipe.md) — the ingredient list.
- [assembly_minimum_information.md](assembly_minimum_information.md) — what each residual term costs.
- `residual_core/docs/limitations.md` · `residual_core/docs/abaqus_validation_roadmap.md` · `STATUS.md`
