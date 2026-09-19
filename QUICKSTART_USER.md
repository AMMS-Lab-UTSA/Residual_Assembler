# Quickstart

This is a one-page guided first run for users. It walks through the three
things Residual_Assembler does, in the order most users need them: sensitivities
of a finished Abaqus analysis, a residual assembled from a model's ingredients,
and sensitivities of a residual you already have. Installation is in
[docs/INSTALL.md](docs/INSTALL.md), every command and option in
[docs/CLI_GUIDE.md](docs/CLI_GUIDE.md), and the graphical interface in
[docs/GUI_GUIDE.md](docs/GUI_GUIDE.md).

Your model stays on your machine: none of these commands makes a network call.

## 1. Sensitivities of a finished Abaqus analysis

You need four files: the input deck `Analysis.inp`, its `Analysis.odb`, the
material compiled as an OTI provider (`OTI_UMAT.obj` with its `Mapping.json`
beside it, built by the companion UMAT-OTI with `umat-oti-provider build`) and
a `sensitivity_request.json`. Then:

```bash
resasm request --model Analysis.inp --odb Analysis.odb \
    --material OTI_UMAT.obj --request sensitivity_request.json --out results
```

A request names the outputs, parameters, region and increments you want:

```json
{
  "outputs": [
    {"name": "loaded_U1", "field": "U", "component": 1, "reduction": "mean"}
  ],
  "parameters": ["E", "SIGY0", "H"],
  "domain": {"nodes": [2, 3, 6, 7]},
  "increments": "LAST"
}
```

Three files are written at the top of `results/`:
`sensitivity_results.json` (values and derivatives), `sensitivity_tables.csv`
(one row per output, increment and parameter) and `run_report.txt` (what was
executed, what was checked, tolerances and limits). The history engine can add
`sensitivity_shares.csv` and, when the request asks for full fields,
`fields.npz`. Everything else (matrices, the linked provider library, the ODB
export) stays in `results/private/`.

`resasm request` handles one bounded J2 model itself and hands every other
readable model (full-size meshes, prescribed displacements, many increments,
any provider) to the history replay engine, `resasm history`, saying so in
`run_report.txt`. Exporting the ODB needs a licensed Abaqus; the sensitivity
computation does not.

**Try it without Abaqus.** The repository ships a small Abaqus J2 beam whose
ODB is already exported to `fields.npz`. From the Residual_Assembler root, with
UMAT-OTI checked out beside it:

```bash
umat-oti-provider build ../UMAT_source_transformation/parameter_sensitivity/models/m3_j2/contract_v2.json \
    --out provider_j2
resasm history --model examples/replay_history/j2_beam/Analysis.inp \
    --fields examples/replay_history/j2_beam/fields.npz \
    --material provider_j2/umat_m3_j2_oti.obj \
    --request examples/replay_history/j2_beam/sensitivity_request.json --out beam_results
```

It prints one summary line and writes the same files, plus
`sensitivity_shares.csv`:

```
history replay: 10 increments, 768 integration points, 4 parameters; yes: max|R_free| = 1.288e-03 N at increment 10 (limit 1.609e-01 N there)
```

Next: the request format and the bounded engine are in
[docs/REQUEST_INTERFACE.md](docs/REQUEST_INTERFACE.md); the history engine, its
mathematics and tolerances in [docs/REPLAY_HISTORY.md](docs/REPLAY_HISTORY.md);
two full-size models in [examples/cantilevers](examples/cantilevers/README.md).
In the GUI this is the **Sensitivity Request** screen.

## 2. Assemble a residual from a model's ingredients

Abaqus never exposes its global residual `R`, but `R` is built from element
residuals, materials, solution fields, loads and constraints. Give those and
Residual_Assembler assembles `R`:

```bash
resasm inspect model.inp                             # what do I have? what is missing?
resasm requirements model.inp --mode stress-driven   # the minimum missing item
resasm init-assembly --model model.inp --solution U.npy
resasm check resasm.yml
resasm run   resasm.yml
```

**`resasm inspect model.inp`** reads the mesh, picks an element backend,
identifies the material and names what it still needs. Real output on a model
shipped in this repository (`sources/permissive/.../HCPnoTwin/Compression111.inp`):

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

**`resasm requirements model.inp --mode stress-driven`** names the one thing to
provide next, never a generic checklist:

```
Minimum missing input:
  provide integration-point stress field S, or an ODB/CSV export (...).
```

**`resasm init-assembly --model model.inp`** writes the recipe and reports what
it inferred, what is still missing and what it can and cannot do. Real output:

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
  constraints            4 *Boundary line(s) read from the mesh
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

Read the **Capability** block. Through this recipe, a C3D8 residual assembles
but is not OTI-differentiated: the element kernels use NumPy float arrays, and
the tool refuses to pretend otherwise. Sensitivities of C3D8 models with a UMAT
come from section 1 instead.

**`resasm check resasm.yml`** repeats that report and stops at the first
blocker. **`resasm run resasm.yml`** does the job.

To assemble `R` directly from an exported stress field, without the recipe:

```bash
abaqus python scripts/extract_odb_fields.py --odb my_job.odb --out fields.json
resasm assemble model.inp --mode stress-driven --fields fields.json --out R.npy
resasm verify   model.inp --fields fields.json     # ||R_free|| ~ 0 at equilibrium
```

`fields.json` is `{"stress_ip": {"<eid>": [[s11,s22,s33,s12,s13,s23], ... per IP]}}`
in Abaqus Voigt order. Other modes: `resasm modes`. Per-mode readiness:
`resasm doctor model.inp`. Backend capabilities: `resasm backends`. The
details are in [docs/residual_assembly_recipe.md](docs/residual_assembly_recipe.md)
and [docs/abaqus_user_path.md](docs/abaqus_user_path.md).

## 3. Sensitivities of a residual you already have

If you can evaluate the global residual `R(u, params)` yourself (a Python
model, a prototype, or a private solver that returns residual coefficients),
start from a template:

```bash
resasm init --template python --out my_case
cd my_case
resasm check resasm.yml
resasm run   resasm.yml
resasm report resasm_output/
```

`init` creates the residual module and its configuration:

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

With `u = 2`, `k = 2`, `f = 16` the run recovers `du/dk = -1/3` and
`du/df = 1/24`, and at order 2 `d²u/dk² = 2/9`.

`--template` is one of `python`, `blackbox`, `blackbox-order2`, `cpp`,
`fortran`. For private code you cannot share, use `blackbox-order2`: your
executable returns Taylor coefficients and the framework only solves the
linear system ([docs/blackbox_order2_contract.md](docs/blackbox_order2_contract.md)).

`check` and `run` accept either dialect: a recipe that names a `mesh:` goes to
the assembly path, a configuration that names a `residual:` goes to the
direct path.

### What you get

```
resasm_output/
  private/   full arrays: residual, tangent, R^(p), U^(p) (= your sensitivities),
             direction_map_order<p>.json, validation_full.json     <- stays local
  public/    summary.md, parameter_ranking.csv, sensitivity_norms.csv,
             validation_summary.json, timing.json                  <- safe to share
```

`public/` holds only norms, rankings, status and timing: no source, no mesh,
no full vectors, no state and no parameter values (names only).

**Coefficients and derivatives (order 2 and above).** OTI produces Taylor
coefficients; the true derivative is `recovery_factor × coefficient` (1 at
order 1, 2 for `d2/dk2`, and so on). Each `.npz` contains both, so read
`U_derivatives`.

### OTILib

The direct Python path needs OTILib, an external GPLv3 library that is not
vendored here: `bash scripts/setup_otilib.sh`, or the no-Conda procedure in
[docs/OTILIB_VENV.md](docs/OTILIB_VENV.md). It does not build natively on
Windows; use WSL (`bash scripts/run_otilib_tests_wsl.sh`). Do **not**
`pip install pyoti`, which is an unrelated package. The black-box path does not
need OTILib. Details:
[residual_core/docs/otilib_integration.md](residual_core/docs/otilib_integration.md).
