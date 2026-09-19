# Command-line guide

Every Residual_Assembler workflow is a subcommand of one program, `resasm`.
This guide describes each subcommand: what it is for, its synopsis, the
options that matter, one worked invocation with its real output (trimmed with
`...`), the files it writes and its exit codes. The last section specifies the
`sensitivity_request.json` format.

All outputs quoted here were measured on 2026-09-18 with the repository's own
examples. Run the commands from the root of the Residual_Assembler checkout,
with the environment of [INSTALL.md](INSTALL.md) active, and send outputs to a
folder outside the repository:

```bash
export WORK="$HOME/resasm_work"; mkdir -p "$WORK"
```

The companion command `umat-oti-provider` (from UMAT-OTI) builds the compiled
material providers that `request`, `history` and `replay` consume; it is
described [at the end](#umat-oti-provider-build).

## Contents

| Group | Subcommand | Purpose |
| --- | --- | --- |
| Sensitivities of a finished analysis | [`request`](#resasm-request) | the four-file interface: `Analysis.inp` + `Analysis.odb` + `OTI_UMAT.obj` + `sensitivity_request.json` |
| | [`history`](#resasm-history) | the history replay engine for any provider; ODB or exported fields; re-equilibration and verification options |
| | [`replay`](#resasm-replay) | the bounded J2 replay of a record or a small model, with an optional finite-difference sweep |
| Inspection | [`inspect`](#resasm-inspect) | what is in a model, which backend, which modes |
| | [`inspect-model`](#resasm-inspect-model) | the ingredients of an assembly recipe: present, inferred, missing |
| | [`requirements`](#resasm-requirements) | what one mode needs and the single missing item |
| | [`doctor`](#resasm-doctor) | inspection plus per-mode readiness; can write a config template |
| | [`modes`](#resasm-modes) | the four assembly modes and their minimum inputs |
| | [`backends`](#resasm-backends) | every registered element backend and its declared limits |
| | [`template`](#resasm-template) | the declared contract of one backend or material |
| Assembly | [`assemble`](#resasm-assemble) | the global residual R (and optionally its tangent) |
| | [`verify`](#resasm-verify) | free residual and reactions of a stress-driven assembly |
| | [`sensitivity`](#resasm-sensitivity) | `T U^(p) = -R^(p)` for a model with parameters, with a finite-difference cross-check |
| Jobs (`resasm.yml`) | [`init`](#resasm-init) | copy a ready-to-run job template, or start the interactive wizard |
| | [`init-assembly`](#resasm-init-assembly) | write an assembly recipe from a model |
| | [`check`](#resasm-check) | is the job ready? One line per requirement |
| | [`run`](#resasm-run) | run the job |
| | [`report`](#resasm-report) | summarise a finished job |

### Conventions

- `resasm --help` lists the subcommands; `resasm <subcommand> --help` shows
  every option.
- `--config FILE` (a `.yml` or `.json` config) is a **global** option. It goes
  before the subcommand: `resasm --config config.json assemble model.json ...`.
- Exit codes follow one pattern: `0` success; `1` the command ran and the
  answer is negative (a check failed, a job is not ready); `2` the command
  could not run with these inputs (missing ingredient, bad path, usage error);
  `3` OTILib was requested and is not installed (`sensitivity` and `run`).
  Each subcommand below lists its measured codes.
- Outputs are split into what you may share and what stays with you. The
  request commands write three public files at the top of `--out` and keep
  full fields, exports and link libraries in `private/`. The job commands
  write `public/` and `private/` folders.

---

## Sensitivities of a finished analysis

### `resasm request`

**Purpose.** The interface for the person who ran the analysis. From the saved
deck, its ODB, the compiled material and a request, it computes the requested
sensitivities without re-running the analysis and without the material
source. Models inside the scope of its bounded single-material J2 engine are
solved there; every other readable model is handed to the history engine
([`history`](#resasm-history)), and the command says so.

**Synopsis.**

```text
resasm request --model MODEL --odb ODB --material MATERIAL --request REQUEST --out OUT
               [--mapping MAPPING] [--abaqus ABAQUS] [--validate]
```

| Option | Meaning |
| --- | --- |
| `--model` | `Analysis.inp`, the deck of the finished analysis |
| `--odb` | `Analysis.odb`; exported with Abaqus Python (`abaqus python`) |
| `--material` | `OTI_UMAT.obj`, the compiled provider |
| `--request` | `sensitivity_request.json` ([format](#the-request-file)) |
| `--out` | a new or empty folder |
| `--mapping` | the completed mapping; by default `<object-stem>.json` or `Mapping.json` beside the object |
| `--abaqus` | the Abaqus launcher used to read the ODB (default `abaqus`) |
| `--validate` | add an independent finite-difference validation with the ORIGINAL routine (history engine: `--verify fd`) |

**Example** ([Example 3](../examples/presentation_request/WALKTHROUGH.md)), in
a folder holding the five input files:

```bash
resasm request --model Analysis.inp --odb Analysis.odb \
    --material OTI_UMAT.obj --request sensitivity_request.json --out results
```

```text
request executed: 4 scalar results; verified=False
```

`verified=False` means that no independent check was requested; the numbers
are computed, not validated. On a model that goes to the history engine
([Example 5](../examples/cantilevers/WALKTHROUGH.md)), the first line of output
names the reason (for example `unsupported *Static options: ['direct']`) and
the last reads `request executed: history engine, 40 increments; verified=False`.

**Outputs.** `sensitivity_results.json`, `sensitivity_tables.csv`,
`run_report.txt` (public); `private/` with the ODB export, the export command
and log, the full-field result and the link library. Runs handed to the
history engine write that engine's outputs, see [`history`](#resasm-history).

**Exit codes.** `0` executed. `2` failed; the message names a category and an
action, `run_report.txt` repeats them, and `private/error_report.txt` holds the
details. Measured categories:

| Category | Example cause |
| --- | --- |
| `output_directory` | `--out` exists and is not empty |
| `odb_export` | Abaqus Python not found (measured with `--abaqus` naming a missing launcher), or the export failed |
| `material_mapping` | no mapping beside the object, or it belongs to another object |
| `request_schema` | the request does not match the format (for example `"parameters": "ALL"` on the bounded engine) |
| `derivative_verification` | `--validate` found a disagreement; on the single-precision ODB of Example 3 its double-precision gate refuses the recorded displacements |

### `resasm history`

**Purpose.** The history replay engine for any provider built by UMAT-OTI:
full-size small-strain C3D8 models, prescribed displacements, many
increments, node and element sets, von Mises outputs, weighted shares and full
fields. It reads the ODB (through Abaqus Python) or an existing export.

**Synopsis.**

```text
resasm history --model MODEL (--odb ODB | --fields FIELDS) --material MATERIAL --request REQUEST --out OUT
               [--mapping MAPPING] [--reequilibrate] [--verify {none,tangent,fd}]
               [--fd-steps FD_STEPS] [--abaqus ABAQUS]
```

| Option | Meaning |
| --- | --- |
| `--odb` / `--fields` | the ODB, or an existing export (`.npz`) made by `residual_core/replay/odb_export_npz.py` (or kept by an earlier run in `private/fields.npz`) |
| `--reequilibrate` | Newton-polish every recorded increment to double-precision equilibrium before taking the sensitivities |
| `--verify tangent` | compare the provider's tangent with central differences of the ORIGINAL UMAT at sampled points |
| `--verify fd` | in addition, re-solve the whole model with the ORIGINAL UMAT at `p (1 +/- h)` and compare every derivative |
| `--fd-steps` | the relative step ladder for `--verify fd` (default `1e-3,3e-4,1e-4,3e-5,1e-5`) |

**Example** ([Example 4](../examples/replay_history/WALKTHROUGH.md)), after
building the J2 provider into `$WORK/provider_j2`:

```bash
resasm history --model examples/replay_history/j2_beam/Analysis.inp \
    --fields examples/replay_history/j2_beam/fields.npz \
    --material "$WORK/provider_j2/umat_m3_j2_oti.obj" \
    --request examples/replay_history/j2_beam/sensitivity_request.json --out "$WORK/beam"
```

```text
history replay: 10 increments, 768 integration points, 4 parameters; yes: max|R_free| = 1.288e-03 N at increment 10 (limit 1.609e-01 N there)
results: .../beam/sensitivity_results.json
```

`run_report.txt` then states each verdict separately:

```text
Status: executed successfully
Command executed: yes: resasm history --model examples/replay_history/j2_beam/Analysis.inp ...
Residual assembled: yes: C3D8 selective-reduced (B-bar), 96 elements, 768 integration points, 585 DOF, 10 increments, sparse assembly
Equilibrium checked: yes: free-DOF residual of the recorded state at all 10 increments
Equilibrium passed: yes: max|R_free| = 1.288e-03 N at increment 10 (limit 1.609e-01 N there)
Tangent available: yes: DDSDDE = dSTRESS/dDSTRAN from the provider's OTI strain directions
Tangent verified: not run
Derivative calculated: yes: total-history du/dp, dRF/dp, dS/dp, dSDV/dp, dMISES/dp for 4 parameters ['E', 'nu', 'SIGY0', 'H'] at 10 increments
Derivative verified: not run
Reference resolved: not applicable (no reference was run)
Abaqus comparison available: yes (primal): replayed S, SDV and RF reproduce the ODB at every integration point and increment; ...
Unsupported feature detected: none
Public and private outputs separated: yes: public ['sensitivity_results.json', 'sensitivity_tables.csv', 'run_report.txt', 'sensitivity_shares.csv']; private private/
```

The last line names every public file the run wrote (here the request asks
for weighted shares, so `sensitivity_shares.csv` is one of them; a
`full_field` request adds `fields.npz`).

With `--reequilibrate --verify fd` (11.6 s here) the lines become
`Tangent verified: yes: max relative error 1.88e-10 vs central FD of the
ORIGINAL UMAT at 36 points ...` and `Derivative verified: yes: whole-model
central FD of the ORIGINAL UMAT re-equilibrated in Python; worst
nonzero-derivative error 1.46e-07 ...`. A derivative counts as verified only
when the finite-difference reference resolved (`Reference resolved: yes`,
adjacent steps agreeing to 1e-4). With coarse steps (`--fd-steps 0.3,0.1`,
measured) the report says `Reference resolved: partially: largest plateau
spread 8.95e-01` and `Derivative verified: not verified: the reference did not
resolve (largest plateau spread 8.95e-01 >= 1e-04); ...`, and
`sensitivity_results.json` records `"verified": false`.

**Outputs.** Public: `sensitivity_results.json` (request, scope, metadata with
parity, tolerances, timings, parameter values and input hashes; every result
with `value`, `derivatives` and `weighted` = `p dQ/dp`),
`sensitivity_tables.csv`, `run_report.txt`, and when requested
`sensitivity_shares.csv` and `fields.npz` (every field and its derivatives at
every increment; it contains the whole solution, so share it only if the model
may be shared). Private: `private/run_details.json`, the ODB export
(`private/fields.npz` when `--odb` was used), its log, and the link library.

**Exit codes.** `0` executed. `2` failed, with the reason on standard error
and in `run_report.txt`; measured examples:
`output directory must be empty or new; refusing to mix results`, and
`Abaqus launcher 'abaqus' not found: exporting Analysis.odb needs Abaqus Python (or pass --fields with an existing export)`.
A deck feature outside the supported subset (several steps, distributed loads,
NLGEOM, other elements, ...) is refused by name
([REPLAY_HISTORY.md](REPLAY_HISTORY.md)).

### `resasm replay`

**Purpose.** The bounded, fingerprint-pinned J2 replay of a replay record
(`resasm_replay_record_v1`), with total-history sensitivities. With `--solve`
it first computes a converged history for a small model with the ORIGINAL
material; with `--verify` it adds whole-path and whole-model finite
differences of the ORIGINAL material. It is the central step of
[Example 6](../examples/bounded_j2_c3d8/WALKTHROUGH.md).

**Synopsis.**

```text
resasm replay RECORD --object OBJECT --contract CONTRACT --out OUT [--solve] [--verify]
```

**Example.**

```bash
resasm replay examples/bounded_j2_c3d8/model.json \
    --object "$WORK/provider_j2/umat_m3_j2_oti.obj" --contract "$WORK/provider_j2/umat_m3_j2_oti.json" \
    --solve --verify --out "$WORK/replay_j2"
```

```text
replay: 7 increments, 8 IPs; total-history du/dp; verified=True
summary: .../replay_j2/public/summary.json
```

**Outputs.** `public/summary.json` (status, parameters, increments, points,
provenance, the verification rows and the largest equilibrium error);
`private/record.json`, `private/result.json`, `private/link/`.

**Exit codes.** `0` completed; `2` failed (`replay failed: <reason>`, and
`public/summary.json` records `"status": "failed"`).

---

## Inspection

### `resasm inspect`

**Purpose.** Read a model (`.inp` or neutral `.json`) and report its elements,
the backend selected for each element type, its materials, the inputs still
missing and the possible assembly modes.

**Synopsis.** `resasm inspect [--detail] model`. `--detail` adds the
per-element-type backend selection with its status and limitations.

**Example.**

```bash
resasm inspect examples/presentation_request/Analysis.inp
```

```text
Model inspection summary
------------------------
Elements:
  - C3D8 elements (x1): supported by solid_c3d8_finite_strain backend

Materials:
  - J2: UMAT adapter required; source file not provided

Required user inputs:
  1. Provide UMAT source or an exported field for material 'J2' (or run stress-driven mode).

Possible modes:
  - formulation (self-contained): not applicable
  - stress-driven residual: available if exported stress/resultant fields are provided (C3D8)
  - material replay: available for C3D8
  - UEL-direct: not applicable
```

When the deck holds keywords that the residual does not apply, `inspect` (and
`doctor`) add a section **Deck keywords present but NOT applied**, one line
per keyword: `*Dsload`, `*Dload`, `*Equation`, `*Amplitude` and `AMPLITUDE=`
references, and any keyword the reader leaves unread. Every command that
reads the deck also prints the same notes as a `DeckKeywordNotApplied`
warning on standard error. Measured on
`tests/cp_c3d8_umat/stress_driven_residual/c3d8_elastic.inp`:

```text
Deck keywords present but NOT applied (the residual omits them):
  - *Elastic: present but not read, so not applied
```

`verify` refuses a deck with any such keyword other than a material
definition (`*Elastic`, `*Plastic`, ...), because the stress-driven residual
takes its stress from the field.

**Outputs.** None (terminal only). **Exit codes.** `0`.

### `resasm inspect-model`

**Purpose.** The assembly-recipe view of a model: which ingredients are
present, which were inferred for you, whether R can be assembled and
OTI-differentiated, and the minimum missing input.

**Synopsis.** `resasm inspect-model [--solution U.npy] [--material FILE] [--param NAME ...] model`

**Example.**

```bash
resasm inspect-model examples/presentation_request/Analysis.inp
```

```text
Residual assembly recipe: Analysis.inp
...
Inferred for you (you did not have to type these):
  constraints            3 *Boundary line(s) read from the mesh
  dof_map                built from the mesh (8 nodes)
  formulation.backend    auto-selected per element type: C3D8 -> solid_c3d8_finite_strain
  ...
  stimuli.loads          1 *Cload line(s) read from the mesh

Capability:
  assemble R           : yes
  OTI-differentiate R  : NO
     blocked: C3D8 -> solid_c3d8_finite_strain

Minimum missing input:
  [missing] fields.solution: the converged solution field U (e.g. U.npy, exported from the ODB)
  [missing] material: a material evaluator (UMAT source, or type + properties) for the solid elements
  [missing] parameters: which parameters to differentiate (e.g. 'MAT.E', 'spring.k')
```

The `OTI-differentiate R: NO` line is deliberate: sensitivities of C3D8 models
with a UMAT go through the compiled provider (`request`, `history`), not
through this recipe.

**Outputs.** None. **Exit codes.** `0` complete; `1` ingredients missing (as here).

### `resasm requirements`

**Purpose.** For one assembly mode, list what is available and name the
single minimum missing input.

**Synopsis.** `resasm requirements --mode MODE [--fields FIELDS] [--subroutine FILE] model`
(`--odb` is another name for `--fields`: a JSON field export, not a binary ODB).

**Example.**

```bash
resasm requirements residual_core/examples/minimal_c3d8_stress_driven/model.json --mode stress-driven
```

```text
Cannot assemble in stress-driven mode.
Available:
  mesh: yes
  formulation backend: yes
  solution field (U / U+rotation / T): yes
  stress / resultant field: no
  one-step deck in scope: yes

Minimum missing input:
  provide integration-point stress field S, or an ODB/CSV export (section resultants N/M/Q for beams/shells; heat flux for thermal).
```

With `--fields residual_core/examples/minimal_c3d8_stress_driven/fields.json`
it answers `Ready to assemble in stress-driven mode.`

A mode is ready only if every element has a backend for it (`formulation
backend`) and, in material replay, a material whose constants are known. A
deck is in scope when it has one step, and `NLGEOM=YES` only with a
finite-strain backend: the general assembly applies all `*Cload` lines of a
deck at once and all its `*Boundary` lines as one list. When the model says
why an input is missing, a `Why:` line follows. Measured:
`--mode formulation` on the cube gives `Why: no registered backend assembles
C3D8 in formulation mode.`; a two-step deck gives `Why: the deck has 2 steps,
and the general assembly would apply the *Cload lines of all of them at once
and their *Boundary lines as one list; resasm history replays a deck step by
step.`

**Outputs.** None. **Exit codes.** `0` in both cases (measured): read the text,
not the code.

### `resasm doctor`

**Purpose.** `inspect` plus a readiness line per mode; optionally writes a
config template to fill in.

**Synopsis.** `resasm doctor [--write-config-template FILE] model`

**Example.**

```bash
resasm doctor residual_core/examples/minimal_c3d8_stress_driven/model.json \
    --write-config-template "$WORK/resasm_config.yml"
```

```text
...
per-mode readiness:
  direct-residual  needs: formulation backend
  formulation      needs: formulation backend
  material-replay  ready
  stress-driven    needs: stress / resultant field

config template written to .../resasm_config.yml
```

Each line is the verdict of `requirements` for that mode, and a mode called
`ready` assembles every element. Here the cube's material replay binds the
section's `E = 210000` and `nu = 0.3` to `solid_c3d8_small_strain`, the C3D8
backend for a small-strain material. The template lists `mode`, `odb`,
`subroutine`, `formulation_policy`,
`material_backend`, `material_parameters` and `options`, each commented.

**Outputs.** The template, if asked for. **Exit codes.** `0`.

### `resasm modes`

**Purpose.** List the assembly modes and their minimum inputs.

```bash
resasm modes
```

```text
assembly modes:
  direct-residual  min inputs: mesh, formulation_backend, element_dofs, uel_routine, deck_scope
  formulation      min inputs: formulation_backend, section_properties, deck_scope
  material-replay  min inputs: mesh, formulation_backend, solution_history, material_model, material_parameters, state_prev, time_increments, deck_scope
  stress-driven    min inputs: mesh, formulation_backend, dof_field, element_field, deck_scope
```

**Exit codes.** `0`.

### `resasm backends`

**Purpose.** An audit of every registered element backend: element types,
DOF types, status, modes, material interface, state, tangent and stated
limitations; and the list of registered materials.

```bash
resasm backends
```

```text
Registered formulation backends
===============================

beam2
  element types : B31, B33, B31_LIKE, BEAM2, FRAME3D
  ...
stress_driven_c3d8
  element types : C3D8
  dof types     : UX, UY, UZ
  status        : verified
  available modes: stress-driven
  material iface : False
  state         : none
  tangent       : none
  limitations   : C3D8 (8-node hex) only; needs exported integration-point stress; no material update; no material tangent (stress-driven)
...
Materials: compressible_neo_hookean, crystal_plasticity, isotropic_elastic, umat
```

Measured backends: `beam2`, `nonlinear_bar1`, `nonlinear_spring1`,
`shell_placeholder` (contract only, not runnable), `solid_c3d8_finite_strain`,
`solid_c3d8_small_strain`, `stress_driven_c3d8`, `truss2`, `uel_direct`
(skeleton). **Exit codes.** `0`.

### `resasm template`

**Purpose.** Print the declared contract of one backend (`--formulation`) or
material (`--material`): required inputs per mode, verification status, tests
and limitations.

```bash
resasm template --formulation stress_driven_c3d8
```

```text
backend: stress_driven_c3d8
  kind:                      formulation
  verification_status:       verified
  supported_element_types:   ['C3D8']
  ...
  required_inputs:           ['coords', 'connectivity', 'dofs', 'stress_ip']
  verification_tests:        ['zero-field', 'patch', 'global-equilibrium']
  ...
```

**Exit codes.** `0`.

---

## Assembly

### `resasm assemble`

**Purpose.** Assemble the global residual `R = F_int - F_ext` of a model in a
given mode and print its size and norms; optionally save it and assemble the
tangent.

**Synopsis.**

```text
resasm assemble [--mode MODE] [--fields FIELDS] [--subroutine FILE] [--tangent] [--out R.npy] model
```

**Example** ([Example 2](../residual_core/examples/minimal_c3d8_stress_driven/WALKTHROUGH.md)):

```bash
M=residual_core/examples/minimal_c3d8_stress_driven
resasm assemble "$M/model.json" --mode stress-driven --fields "$M/fields.json" --out "$WORK/R_cube.npy"
```

```text
assembled mode 'stress-driven': ndof=24  ||R||=7.071068e+01  max|R|=2.500000e+01
residual saved to .../R_cube.npy
```

A finite-strain model with its config
([Example 7](../examples/finite_strain_c3d8/WALKTHROUGH.md)):
`resasm --config config.json assemble model.json --mode material-replay --tangent`
printed `ndof=36  ||R||=1.744643e+00  max|R|=5.938249e-01`.

**Outputs.** `--out` saves R as a NumPy `.npy` array (node by node, DOF by DOF).

**Exit codes.** `0` assembled, with every element of the model. `2` the mode
is not runnable (the requirements report is printed instead of a residual,
with its `Why:` line; measured on `minimal_truss --mode material-replay`:
`Why: no registered backend assembles T3D2 in material-replay mode.`), or the
`--fields` file cannot be used. One line then names the file and the reason, measured:
`ERROR: cannot read the field export .../field.json: No such file or directory`,
and for a file in another layout (here the model file)
`ERROR: field export .../model.json: key 'schema' is not an element id; expected {"stress_ip": ...}`.

### `resasm verify`

**Purpose.** For a stress-driven assembly, split R into the free part (which
should vanish at equilibrium) and the reactions at the constrained DOFs, and
test the free part against an absolute tolerance.

**Synopsis.** `resasm verify [--fields FIELDS] [--atol ATOL] model` (default
`--atol 1e-6`, in model force units).

**Example.**

```bash
resasm verify "$M/model.json" --fields "$M/fields.json"
```

```text
stress-driven verification:
  ndof                 = 24
  ||R_free||           = 7.071068e+01   (should be ~0 at equilibrium)
  ||reaction (at BC)|| = 0.000000e+00
  absolute tolerance   = 1.000000e-06 (model force units)
  equilibrium          = FAIL
  reaction reference   = NOT CHECKED (assembled reactions only)
```

`FAIL` is correct here: the example cube carries a stress but no supports or
loads. Reactions are printed but not compared with a reference.

**Exit codes.** `0` equilibrium within `--atol`; `1` equilibrium failed (as
here); `2` verification could not run: a `--fields` file that is missing or
in another layout (the same one-line `ERROR:` as `assemble`), a mode that is
not runnable (for example a deck with several steps or `NLGEOM=YES`), or a
deck with loads or constraints the residual does not apply. The last prints
`Cannot verify equilibrium: the deck asks for what the residual does not
apply, so its equilibrium cannot be verified:` followed by the keywords
(measured with `*Dsload`, `*Dload`, `*Amplitude` and `AMPLITUDE=` added to the
one-element deck).

### `resasm sensitivity`

**Purpose.** For a model whose backend depends on named parameters, build the
hypercomplex right-hand sides `R^(p)` and solve `T U^(p) = -R^(p)` order by
order, with a finite-difference cross-check of the first order.

**Synopsis.**

```text
resasm sensitivity [--params FILE] [--param NAME ...] [--mode MODE] [--order N]
                   [--backend otilib|dual1] [--no-fd] [--out PREFIX] model
```

| Option | Meaning |
| --- | --- |
| `--params` | a JSON or YAML file with `parameters`, `mode`, `order`, `backend` and optional `expected` values |
| `--param` | a parameter `group.key` (repeatable); overrides the file |
| `--order` | derivative order (OTILib supports 1 and higher) |
| `--backend` | `otilib` (default, production) or `dual1` (a legacy first-order smoke test for backends whose parameters are section entries; it refuses the finite-strain formulation) |
| `--no-fd` | skip the finite-difference cross-check |
| `--out` | save the package as `PREFIX.{json,npz,md}` plus `PREFIX_residual.*`, `PREFIX_sensitivity.*`, `PREFIX_state.*`, `PREFIX_validation.*` |

**Example** (the built-in cubic spring `R = k u^3 - f`, needs OTILib):

```bash
S=residual_core/examples/minimal_nonlinear_spring_sensitivity
resasm sensitivity "$S/model.json" --param spring.k --param spring.f --order 2 --out "$WORK/spring_sens"
```

```text
hypercomplex backend     : otilib
basis count (m)          : 2
truncation order (nt)    : 2
total coefficients (N)   : 6
sensitivity analysis (mode=formulation, order=2)
  residual norm ||R_free|| = 0.000000e+00
  tangent source           = backend-assembled
  parameters               = ['spring.k', 'spring.f']
  R^(1) shape              = (1, 2)  (2 directions)
  R^(2) shape              = (1, 3)  (3 directions)
  hypercomplex ready       = True
  solved derivative orders : [1, 2]
  order 1:
    d^1/e1         = -3.333333e-01   [FD -3.333330e-01, rel 1.00e-06]
    d^1/e2         = +4.166667e-02   [FD +4.166667e-02, rel 1.40e-10]
  order 2:
    d^2/e1^2       = +2.222222e-01
    d^2/e1*e2      = -6.944444e-03
    d^2/e2^2       = -1.736111e-03
  saved package -> .../spring_sens.{json,npz,md}
```

The printed values are recovered derivatives (`e1^2` includes the factor 2);
for a model with several DOFs they are the Euclidean norms of the solution
sensitivity columns. With `--params "$S/params.json"`, which lists an
`expected` value, the line also shows `[analytic -3.333333e-01, rel 0.00e+00]`.
The `[FD ...]` column compares every solved first-order column with finite
differences of re-solved equilibria; `rel` is the relative error of the whole
column. If any `rel` reaches `1e-4`, the command prints
`finite-difference check FAILED for <labels>: relative error >= 0.0001; these
derivatives are not confirmed` and exits 1. A backend that cannot produce a
derivative refuses instead of printing a zero, exit 2. Measured:
`--param spring.zzz` gives `ERROR: backend 'otilib' cannot differentiate with
respect to 'spring.zzz': material 'spring' has no section entry 'zzz'
(entries: f, k)`, and `--backend dual1` on the finite-strain model of
[Example 7](../examples/finite_strain_c3d8/WALKTHROUGH.md) gives
`ERROR: backend='dual1' cannot differentiate the finite-strain C3D8
formulation (solid_c3d8_finite_strain): only the OTILib backend seeds its
material constants. Use backend='otilib'.`

**Exit codes.** `0` solved, and every first-order derivative agrees with the
finite differences (unless `--no-fd`); `1` a finite-difference check failed;
`2` the mode is not runnable (measured: `--mode stress-driven` without a field
export prints `Cannot assemble in stress-driven mode.` and the missing input),
the system is not ready, or the backend cannot run, with the cause on the
`ERROR:` line (measured on the stress-driven cube, its field attached with a
config `odb:` entry: `ERROR: cannot solve R(u) = 0 in stress-driven mode: no
element assembles a tangent dR/du there, so there is no Newton step`); `3`
OTILib was requested and is not installed, with an installation hint.

---

## Jobs: `resasm.yml`

A job is a folder with a `resasm.yml` and whatever it refers to. Paths in
`resasm.yml` are resolved against the folder that holds it, so the commands
work from any directory. `check` and `run` accept both dialects of
`resasm.yml`: a file with `mesh:` is an assembly recipe, a file with
`residual:` is a direct or black-box residual job
([simple_config_contract.md](simple_config_contract.md)).

### `resasm init`

**Purpose.** Start a job: copy a ready-to-run template into a new folder, or,
without `--template`, run an interactive wizard in the terminal that writes a
`resasm.yml`.

**Synopsis.** `resasm init [--template {python,blackbox,blackbox-order2,cpp,fortran}] [--out OUT] [--force]`

| Template | What you provide |
| --- | --- |
| `python` | `residual(u, params)` in Python; derivatives by OTILib at any order |
| `blackbox` | an executable that returns first-order residual coefficients |
| `blackbox-order2` | an executable that returns residual coefficients up to order 2; no OTILib needed on this side |
| `cpp`, `fortran` | the residual in C++ or Fortran, with a build file |

**Example.**

```bash
resasm init --template python --out "$WORK/spring"
```

```text
created .../spring (from template 'python')
  README.md
  resasm.yml
  solution.npy
  user_residual.py

Next:
  cd .../spring
  resasm check resasm.yml
  resasm run resasm.yml
```

**Exit codes.** `0` created; `2` the folder exists
(`ERROR: ... already exists (use --force to overwrite)`).

### `resasm init-assembly`

**Purpose.** Write an assembly recipe (`resasm.yml` with `mesh:`) from a
model, and report what was inferred and what is still missing.

**Synopsis.** `resasm init-assembly [--model MODEL] [--solution U.npy] [--material FILE] [--param NAME ...] [--name NAME] [--out resasm.yml]`

**Example.**

```bash
resasm init-assembly --model "$PWD/examples/presentation_request/Analysis.inp" --name one_element \
    --out "$WORK/one_element_recipe.yml"
```

It prints `wrote .../one_element_recipe.yml`, then the same report as
[`inspect-model`](#resasm-inspect-model), and `Next: resasm check ...`. Give
the model as an absolute path, as here: the recipe records the path as typed,
and `check` and `run` may read it from another folder.

**Exit codes.** `0` written (even when ingredients are still missing).

### `resasm check`

**Purpose.** Check that a job is ready to run: one `[ok]`, `[warn]` or
`[fail]` line per requirement, stopping at the first actionable failure.

**Example** (the `python` template):

```bash
resasm check "$WORK/spring/resasm.yml"
```

```text
[ok] config parsed: spring_demo (order 2, backend otilib)
[ok] loaded solution vector: shape (1,)
[ok] loaded parameter map: 2 parameters (k, f)
[ok] residual evaluated: shape (1,)
[ok] tangent loaded (python): shape (1, 1)
[ok] OTILib available: order 2, basis 2
[ok] residual norm on free DOFs is 0.000e+00
[ok] RHS order 1 generated: shape (1, 2)
[ok] sensitivity solve completed
```

Without OTILib the sixth line reads
`[fail] OTILib backend requested but genuine OTILib was not found.`
For an incomplete assembly recipe `check` prints the recipe report with its
`[missing]` items. For a file that does not exist it prints
`[fail] Config file not found: ...` and a minimal `resasm.yml` to start from.

**Exit codes.** `0` ready; `1` not ready.

### `resasm run`

**Purpose.** Run the job: evaluate the residual, generate the right-hand
sides, solve every order, run the configured validation and write the output
package.

**Example.**

```bash
resasm run "$WORK/spring/resasm.yml"
```

```text
run complete: {'parameters': ['k', 'f'], 'order': 2, 'tangent_source': 'python', 'residual_free_norm': 0.0, 'orders_solved': [1, 2]}
  private outputs -> .../spring/resasm_output/private
  public  outputs -> .../spring/resasm_output/public
  read: .../spring/resasm_output/public/summary.md
```

The `blackbox-order2` template, run the same way, returned
`du/dk = -0.66666667`, `du/df = 0.02083333`, `d2u/dk2 = 0.555555556`,
`d2u/dk df = -0.00694444444` and `d2u/df2 = -0.000434027778` (read from
`private/solution_sensitivities_order{1,2}.npz`), which are -2/3, 1/48, 5/9,
-1/144 and -1/2304, the closed form for `R = k^2 u^3 - f` at `k = 2`,
`f = 32`, `u = 2`.

**Outputs.** `resasm_output/` next to `resasm.yml` (or `output: dir:` in the
file): `public/` (`summary.md`, `sensitivity_norms.csv`,
`parameter_ranking.csv`, `validation_summary.json`, `timing.json`) and
`private/` (solution sensitivities of every order with coefficients, recovered
derivatives and recovery factors; right-hand sides; residual; tangent;
direction maps; metadata).

**Exit codes.** `0` complete. `2` an incomplete assembly recipe
(`ERROR: the assembly recipe is not complete:` followed by the report).
`3` the job needs OTILib (the `python` template and assembly recipes) and it
is not installed: `ERROR: OTILib backend requested but genuine OTILib was not
found.`, followed by how to install it.

### `resasm report`

**Purpose.** Summarise a finished job from its output folder, or from its
`resasm.yml`: then it reads the folder that job writes to (`resasm_output/`
beside the file, or `output: dir:`), as the GUI's **Read report** does.

```bash
resasm report "$WORK/spring/resasm_output"
resasm report "$WORK/spring/resasm.yml"      # the same report
```

```text
Sensitivity run report: .../spring/resasm_output
  residual norm (free) : 0.0
  tangent source       : python
  parameters           : k, f
  derivative orders    : [1, 2]
  validation status    : ok
  private outputs      : .../spring/resasm_output/private
  public  outputs      : .../spring/resasm_output/public
```

**Exit codes.** `0`; `2` if the folder has no `public/validation_summary.json`
(`ERROR: Cannot read report file .../public/validation_summary.json: ...`).

---

## `umat-oti-provider build`

Part of UMAT-OTI. It compiles the ORIGINAL UMAT and its OTI version into one
object and writes the completed mapping beside it:

```bash
umat-oti-provider build ../UMAT_source_transformation/parameter_sensitivity/models/m3_j2/contract_v2.json \
    --out "$WORK/provider_j2"
```

It prints a JSON object with the paths of `umat_m3_j2_oti.obj`,
`umat_m3_j2_oti.json` and the build folder (measured 5.5 s for J2, 8.6 s for
the FCC crystal model). Hand over only the `.obj` and the `.json`; they may be
renamed to `OTI_UMAT.obj` and `Mapping.json`, since the mapping records the
object's SHA-256. The build folder contains generated sources and stays with
the material owner. Options: `--compiler`, `--regular-object NAME.obj` (also
publish the ORIGINAL UMAT compiled unchanged), `--abaqus-toolchain`,
`--abaqus`.

---

## The request file

`sensitivity_request.json` says which quantities to differentiate, with
respect to which parameters, where, and at which increments. Both engines
read the same core format; the history engine accepts a few extensions.
`resasm request` keeps a request on its bounded engine only when the model is
inside that engine's scope **and** the request uses the core format; anything
else goes to the history engine.

### Keys

| Key | Required | Value |
| --- | --- | --- |
| `outputs` | yes | a non-empty list of outputs (below) |
| `parameters` | yes | a list of unique parameter names from the mapping, e.g. `["E", "SIGY0", "H"]`; the history engine also accepts `"ALL"` (the provider's order) |
| `domain` | yes | the default region for outputs without their own `domain`, e.g. `{"nodes": "ALL", "elements": "ALL"}` |
| `increments` | yes | `"LAST"`, `"ALL"` or a list of unique one-based increment numbers, e.g. `[2, 4]`. The whole history before them is always replayed |
| `weighted_shares` | no (history) | `{"field": "MISES", "domain": {"elements": "ALL"}}`; `field` may be `MISES` (default), `S` or `SDV`. Writes `sensitivity_shares.csv` |
| `full_field` | no (history) | `true` writes `fields.npz` with every field and its derivatives at every increment |

### One output

```json
{"name": "tip_RF2", "field": "RF", "component": 2, "reduction": "sum", "domain": {"nset": "TIP"}}
```

| Key | Value |
| --- | --- |
| `name` | a unique non-empty string; it labels the rows of the tables |
| `field` | `U` (displacement) or `RF` (reaction) at nodes; `S` (stress), `SDV` (state variables) or `MISES` (von Mises stress, history engine) at integration points |
| `component` | one-based, or `"ALL"`. `U`, `RF`: 1 to 3. `S`: 1 to 6 in the order 11, 22, 33, 12, 13, 23. `SDV`: 1 to the number of state variables (the bounded engine: `SDV1`, the equivalent plastic strain). `MISES`: 1 |
| `reduction` | how the selected values become one number (below) |
| `domain` | optional; overrides the request's `domain` |

### Reductions

| Reduction | Result | Engines |
| --- | --- | --- |
| `component` | the value at exactly one location (one node or one point, one component) | both |
| `sum` | the plain sum (not a spatial integral) | both |
| `mean` | the plain mean (not volume-weighted) | both |
| `volume_mean` | the mean weighted by the integration-point volumes `det(J) w_q`; integration-point fields only | history |
| `L2` | the Euclidean norm (not von Mises); a zero field is refused as not differentiable | both |
| `max`, `min` | the signed extreme. If it is attained at several locations whose derivatives differ, the output is refused as not differentiable | `max`: both; `min`: history |

### Domains

| Key | Meaning | Engines |
| --- | --- | --- |
| `nodes` | `"ALL"` or a list of node ids (for `U`, `RF`) | both |
| `elements` | `"ALL"` or a list of element ids (for `S`, `SDV`, `MISES`); all eight integration points unless `points` is given | both |
| `nset`, `elset` | a node or element set of the `.inp` | history |
| `points` | a list of integration points, 1 to 8 | history |

### Examples

The core format, as in `examples/presentation_request/sensitivity_request.json`
(four outputs at the last increment; all four J2 parameters):

```json
{
  "outputs": [
    {"name": "mean_loaded_U1", "field": "U", "component": 1, "reduction": "mean", "domain": {"nodes": [2, 3, 6, 7]}},
    {"name": "support_RF1", "field": "RF", "component": 1, "reduction": "sum", "domain": {"nodes": [1, 4, 5, 8]}},
    {"name": "mean_S11", "field": "S", "component": 1, "reduction": "mean"},
    {"name": "mean_EQPLAS", "field": "SDV", "component": 1, "reduction": "mean"}
  ],
  "parameters": ["E", "nu", "SIGY0", "H"],
  "domain": {"nodes": "ALL", "elements": "ALL"},
  "increments": "LAST"
}
```

The history extensions, as in `examples/cantilevers/j2_request.json`:

```json
{
  "outputs": [
    {"name": "tip_RF2", "field": "RF", "component": 2, "reduction": "sum", "domain": {"nset": "TIP"}},
    {"name": "tiptop_U1", "field": "U", "component": 1, "reduction": "component", "domain": {"nodes": [833]}},
    {"name": "mises_mean", "field": "MISES", "component": 1, "reduction": "volume_mean", "domain": {"elements": "ALL"}},
    {"name": "mises_root_max", "field": "MISES", "component": 1, "reduction": "max", "domain": {"elset": "ROOTEL"}},
    {"name": "S11_e1_ip1", "field": "S", "component": 1, "reduction": "component", "domain": {"elements": [1], "points": [1]}},
    {"name": "eqplas_max", "field": "SDV", "component": 1, "reduction": "max", "domain": {"elements": "ALL"}}
  ],
  "parameters": "ALL",
  "domain": {"nodes": "ALL", "elements": "ALL"},
  "increments": "ALL",
  "weighted_shares": {"field": "MISES", "domain": {"elements": "ALL"}},
  "full_field": true
}
```

### Rules worth knowing

- Unknown keys, fields, reductions or ids are refused; nothing is silently
  dropped.
- `"parameters": "ALL"` is a history-engine extension. On a model that stays
  on the bounded engine (Example 3) it was refused with `Category:
  request_schema` (measured); give the list of names instead.
- Adding a `MISES` output (or any other extension) to a request moves it to
  the history engine. Measured on the one-element model of Example 3: a
  `volume_mean` of `MISES` was accepted there and returned 299.99999955803224
  with derivatives of order 1e-14 or smaller, as the closed form requires.
- The parameter value used for `weighted` (`p dQ/dp`) and for the shares is
  the constant in the `.inp` at the PROPS index that the mapping gives.
- Only first derivatives with respect to material parameters are computed.
  Load, boundary-condition and shape sensitivities are not supported.

For the mathematics, tolerances and the supported deck subset, see
[REPLAY_HISTORY.md](REPLAY_HISTORY.md) and
[REQUEST_INTERFACE.md](REQUEST_INTERFACE.md). The worked examples are in
[examples/README.md](../examples/README.md); the GUI equivalents are in
[GUI_GUIDE.md](GUI_GUIDE.md).
