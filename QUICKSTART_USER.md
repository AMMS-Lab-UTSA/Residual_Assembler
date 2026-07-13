# Quickstart (users)

You should not have to provide the global residual R. Abaqus never exposes it.
Provide the ingredients — mesh, material, solution field, loads, BCs — and
Residual_Assembler assembles R for you. Your model stays on your machine.

**Install:** `pip install -e .` (exposes `resasm`). Assembly needs nothing else.

## The assembly path (start here)

```
resasm inspect model.inp                             # what do I have? what's missing?
resasm requirements model.inp --mode stress-driven   # the minimum missing item
resasm init-assembly --model model.inp --solution U.npy
resasm check resasm.yml
resasm run   resasm.yml
```

**`resasm inspect model.inp`** reads the mesh, picks an element backend, identifies
the material, and names what it still needs. Real output on a model shipped in this
repo (`sources/permissive/.../HCPnoTwin/Compression111.inp`):

```
Model inspection summary
------------------------
Elements:
  - C3D8 elements (x125): supported by solid_c3d8_finite_strain backend

Materials:
  - CPuranium: UMAT adapter required; source file not provided

Required user inputs:
  1. Provide UMAT source or an exported field for material 'CPuranium' (or run stress-driven mode).

Possible modes:
  - formulation (self-contained): not applicable
  - stress-driven residual: available if exported stress/resultant fields are provided (C3D8)
  - material replay: available for C3D8
  - UEL-direct: not applicable
```

**`resasm requirements model.inp --mode stress-driven`** names the *one* thing to
provide next — never a generic checklist:

```
Minimum missing input:
  provide integration-point stress field S, or an ODB/CSV export (...).
```

**`resasm init-assembly --model model.inp`** writes the recipe and immediately tells
you what it inferred, what is still missing, and — honestly — what it can and cannot
do. Real output:

```
Residual assembly recipe: assembly_job
----------------------------------------------
Mesh:
  file            : .../HCPnoTwin/Compression111.inp
  nodes/elements  : 216 / 125
  element types   : C3D8
Ingredients:
  solution field  : -- MISSING --
  material        : -- MISSING --
  element fields  : (not given)
  parameters      : -- MISSING --

Inferred for you (you did not have to type these):
  constraints            4 *Boundary block(s) read from the mesh
  dof_map                built from the mesh (216 nodes)
  formulation.backend    auto-selected per element type: C3D8 -> solid_c3d8_finite_strain
  mesh.element_types     from Compression111.inp (C3D8 x125)
  sensitivity.backend    default otilib
  sensitivity.order      default 1

Capability:
  assemble R           : yes
  OTI-differentiate R  : NO
     blocked: C3D8 -> solid_c3d8_finite_strain

Minimum missing input:
  [missing] fields.solution: the converged solution field U (e.g. U.npy, exported from the ODB)
  [missing] material: a material evaluator (UMAT source, or type + properties) for the solid elements
  [missing] parameters: which parameters to differentiate (e.g. 'MAT.E', 'spring.k')
```

Read the **Capability** block. For C3D8 today: **R assembles, but R cannot be
OTI-differentiated** — the element kernels use numpy float arrays. The tool refuses
to pretend otherwise. See [README.md](README.md#status).

**`resasm check resasm.yml`** re-runs that report and stops at the first blocker.
**`resasm run resasm.yml`** does the job.

To assemble R directly from an exported stress field, without the recipe:

```
abaqus python scripts/extract_odb_fields.py --odb my_job.odb --out fields.json
resasm assemble model.inp --mode stress-driven --fields fields.json --out R.npy
resasm verify   model.inp --fields fields.json     # ||R_free|| ~ 0 at equilibrium
```

`fields.json` is `{"stress_ip": {"<eid>": [[s11,s22,s33,s12,s13,s23], ... per IP]}}`
in Abaqus Voigt order. Other modes: `resasm modes`. Per-mode readiness:
`resasm doctor model.inp`. Backend capabilities: `resasm backends`.

## If you already have `R(u, params)`

Secondary path — the shortcut. Only usable if you can evaluate the *global* residual
yourself (an Abaqus user cannot). It is also the smallest end-to-end example, and it
is a working sensitivity route.

```
resasm init --template python --out my_case
cd my_case
resasm check resasm.yml
resasm run   resasm.yml
resasm report resasm_output/
```

`init` creates two files:

```python
# user_residual.py
def residual(u, params, state=None, time=None):
    k = params["k"]; f = params["f"]
    return [k * u[0] ** 3 - f]        # R(u, params)

def tangent(u, params, state=None, time=None):
    return [[3.0 * params["k"] * u[0] ** 2]]   # dR/du
```

```yaml
# resasm.yml
problem:     { name: spring_demo, unknowns: 1 }
residual:    { type: python, module: user_residual.py, function: residual }
tangent:     { type: python, function: tangent }
parameters:  { k: 2.0, f: 16.0 }      # the design variables
solution:    { file: solution.npy }   # your converged u
sensitivity: { order: 2, backend: otilib }
validation:  { rhs_finite_difference_check: true }
```

`--template` is one of `python`, `blackbox`, `blackbox-order2`, `cpp`, `fortran`.
For private code you cannot share, use `blackbox-order2`
([docs/blackbox_order2_contract.md](docs/blackbox_order2_contract.md)) — that is the
working sensitivity route today.

`check` and `run` accept **either** dialect: a recipe that names a `mesh:` routes to
the assembly path, a config that hands over a `residual:` routes to the direct path.

## What you get

```
resasm_output/
  private/   full arrays: residual, tangent, R^(p), U^(p) (= your sensitivities),
             direction_map_order<p>.json, validation_full.json     <- stays local
  public/    summary.md, parameter_ranking.csv, sensitivity_norms.csv,
             validation_summary.json, timing.json                  <- safe to share
```

`public/` holds only norms, rankings, status and timing — no source, no mesh, no full
vectors, no state, no parameter *values* (names only).

**Coefficients vs derivatives (order ≥ 2):** OTI gives Taylor *coefficients*; the true
derivative is `recovery_factor × coefficient` (1 at order 1, 2 for `d2/dk2`, …). Each
`.npz` ships both — **just read `U_derivatives`.**

## OTILib

The sensitivity paths need OTILib (GPLv3, external, not vendored):
`bash scripts/setup_otilib.sh`. It does not build natively on Windows — use WSL
(`bash scripts/run_otilib_tests_wsl.sh`). Do **not** `pip install pyoti` (unrelated
package). Details:
[residual_core/docs/otilib_integration.md](residual_core/docs/otilib_integration.md).
