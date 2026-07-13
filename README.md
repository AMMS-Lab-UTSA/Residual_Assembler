# Residual_Assembler

**You should not have to provide R. You provide the ingredients, and
Residual_Assembler builds R.**

Abaqus hides the global residual. It owns it internally and never exposes it — so
"just write your residual function" is not something an Abaqus user can do. But R
is not magic. It is assembled from element residuals, materials, solution fields,
loads and constraints. So those are what we ask for.

```
R(u, a) = F_internal(u, a, q) - F_external(a, t) + F_constraints(u, t)

R_e     = ∫_Ωe B^T σ(u, a, q) dΩ - f_e^ext
```

> **One model recipe. One converged solution. One parameter list. Sensitivities out.**

Your model never leaves your machine.

**Read the Status section before you plan around this.** Today the assembly beats
are verified for C3D8: you can assemble and verify R from exported Abaqus
ingredients. The *sensitivity* beat does **not** yet run through the C3D8 assembly
path — see [Status](#status) for exactly what works and what does not.

## What you provide (the ingredients)

| Ingredient | Where it comes from |
|---|---|
| **Mesh** — nodes, connectivity, element types | your `.inp` (parsed for you) |
| **Element formulation** | auto-selected from the element type (`resasm inspect`) |
| **Material** | your `.inp` material, a UMAT source, or a built-in |
| **Solution field** — converged `u` (U, U+rotation, T) | your solve / ODB export |
| **State / history** — STATEV, internal variables | ODB export (material-replay) |
| **Loads, BCs, stimuli** | your `.inp` (`*Cload`, `*Boundary`) |
| **Parameters** — the design variables | you name them |
| **Tangent source** | assembled, exported, or your own |

You do not hand over all of these at once. The tool tells you the **single minimum
missing item** for the mode you asked for, and nothing more.

## What we build

The global residual `R`, assembled element by element, plus the free/prescribed
split and the reactions at constrained DOFs. `R = F_internal - F_external`
(`residual_core/core/assembler.py`).

## The Abaqus flow

Start with nothing but a `.inp`. This is a real run against a model shipped in this
repo:

```
resasm inspect sources/permissive/ngrilli_Oxford_Crystal_Plasticity/ExampleInputFiles/HCPnoTwin/Compression111.inp
```

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

It read the mesh, found the elements, picked a backend, identified the material,
and told you the one thing it still needs. Now ask what a specific mode requires:

```
resasm requirements <model.inp> --mode stress-driven
```

```
Cannot assemble in stress-driven mode.
Available:
  mesh: yes
  solution field (U / U+rotation / T): yes
  stress / resultant field: no

Minimum missing input:
  provide integration-point stress field S, or an ODB/CSV export (section resultants N/M/Q for beams/shells; heat flux for thermal).
```

One missing item, named. Then build the recipe:

```
resasm init-assembly --model <model.inp> --solution U.npy
resasm check resasm.yml
resasm run   resasm.yml
```

`init-assembly` writes the recipe and immediately reports what it inferred, what is
still missing, and — honestly — what it can and cannot do. Real output:

```
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

The **Capability** block is the honest part: for C3D8 today R assembles, but R cannot
be OTI-differentiated. The tool refuses to pretend otherwise. See [Status](#status).

To assemble R directly from an exported stress field, without the recipe:

```
abaqus python scripts/extract_odb_fields.py --odb my_job.odb --out fields.json

resasm assemble <model.inp> --mode stress-driven --fields fields.json --out R.npy
resasm verify   <model.inp> --fields fields.json
```

`assemble` prints `ndof`, `||R||` and `max|R|`. `verify` splits the residual and
reports `||R_free||` (should be ~0 at equilibrium) against the reaction norm at the
constrained DOFs.

The field export is JSON:
`{"stress_ip": {"<eid>": [[s11,s22,s33,s12,s13,s23], ... per IP]}}` — Abaqus Voigt
order. `resasm assemble` reads `stress_ip`; the same file also carries `U`, `RF`
and `SDV` for the other paths.

Other commands that exist today: `resasm inspect-model <model>` (the recipe-shaped
inspector), `resasm doctor <model>` (per-mode readiness, and
`--write-config-template` to emit a config), `resasm backends` (what every backend
declares), `resasm modes`, `resasm template --formulation <name>`.

## The three paths

### Path A — assembly from ingredients (the main path)

You give the mesh, material, solution field, state, loads and BCs. We assemble R.
This is the path for an Abaqus-like solver that never exposes its residual.

Modes (`resasm modes`), each with its own minimum inputs:

| Mode | Minimum inputs |
|---|---|
| `stress-driven` | mesh, DOF solution field, element field (IP stress / resultants) |
| `material-replay` | mesh, solution history, material model, PROPS, previous STATEV, time increments |
| `direct-residual` | mesh, element DOF layout, a UEL-like callable |
| `formulation` | a registered formulation backend, section properties |

Recipe and the exact minimum information:
[docs/residual_assembly_recipe.md](docs/residual_assembly_recipe.md),
[docs/assembly_minimum_information.md](docs/assembly_minimum_information.md),
[docs/abaqus_user_path.md](docs/abaqus_user_path.md).

### Path B — black-box

Your private executable reads a `request.json` and writes back residual
coefficients (order ≥ 2). We never see your code, mesh, or material model; we only
solve `T U^(p) = -R^(p)`.

```
resasm init --template blackbox-order2 --out my_case
```

Contract: [docs/blackbox_order2_contract.md](docs/blackbox_order2_contract.md).
This is a working sensitivity route today.

### Path C — smallest toy example (if you already have R)

**This is the shortcut, not the collaborator workflow.** It only applies if you can
already evaluate the *global* residual yourself — which an Abaqus user cannot. It
exists because it is the smallest thing that demonstrates the sensitivity machinery
end to end, and because some users genuinely do own their residual.

You write `residual(u, params)`:

```python
def residual(u, params, state=None, time=None):
    k = params["k"]
    f = params["f"]
    return [k * u[0] ** 3 - f]        # cubic spring:  R = k u^3 - f


def tangent(u, params, state=None, time=None):
    k = params["k"]
    return [[3.0 * k * u[0] ** 2]]    # dR/du
```

Use ordinary arithmetic. The framework calls `residual()` once with plain floats and
once with hypercomplex (OTI) numbers, and reads the derivatives off the result.

```
resasm init --template python --out my_case
cd my_case
resasm check resasm.yml
resasm run   resasm.yml
resasm report resasm_output/
```

With `u = 2`, `k = 2`, `f = 16` this recovers `du/dk = -1/3` and `du/df = 1/24`.

`--template` is one of `python`, `blackbox`, `blackbox-order2`, `cpp`, `fortran`.
`resasm check` is a readiness report that stops at the first actionable `[fail]`.

`check` and `run` accept **either** dialect of `resasm.yml`: a recipe that names a
`mesh:` routes to the assembly path (A), a config that hands over a `residual:`
routes to the direct path (B / C).

## The sensitivity contract

Where a path does produce sensitivities, the object is

```
T U^(p) = -R^(p)
```

| symbol | meaning |
|---|---|
| `T` | tangent `dR/du` at the converged real solution |
| `R^(p)` | p-th order residual coefficients from the OTI evaluation |
| `U^(p)` | p-th order solution sensitivity coefficients — what you wanted |

> **Coefficients vs derivatives.** OTI produces Taylor *coefficients*. The true
> partial derivative is `recovery_factor × coefficient`, with
> `recovery_factor = Π_i (κ_i!)` — `1` at order 1, `2!` for `d2/dk2`, and so on. You
> do not have to apply it: every order exports **both** (`U_coefficients` /
> `U_derivatives`) plus `direction_map_order<p>.json`. Public reports quote
> recovered derivatives.

Outputs split into `private/` (full arrays — stays on your machine) and `public/`
(norms, rankings, status, timing — the folder you can share). See
[docs/privacy_contract.md](docs/privacy_contract.md) and
[docs/output_objects.md](docs/output_objects.md).

## Install

```
pip install -e .            # exposes the `resasm` command
```

Requires Python ≥ 3.9 and numpy. `PyYAML` is optional (`pip install -e .[yaml]`); a
minimal YAML reader is bundled.

Path A (assembly, inspection, verification) needs nothing else.

**OTILib** is the hypercomplex engine behind the sensitivity paths. It is an
**external GPLv3 dependency, not vendored here**:

```
bash scripts/setup_otilib.sh          # builds https://github.com/mauriaristi/otilib.git
```

It does not build natively on Windows — use WSL
(`bash scripts/run_otilib_tests_wsl.sh`). Do **not** `pip install pyoti`; that PyPI
name is an unrelated package. See
[residual_core/docs/otilib_integration.md](residual_core/docs/otilib_integration.md).
The black-box path does not require OTILib on our side.

## Docs

Assembly (Path A):
- [docs/residual_assembly_recipe.md](docs/residual_assembly_recipe.md) — how R is assembled from ingredients
- [docs/assembly_minimum_information.md](docs/assembly_minimum_information.md) — the minimum information per mode
- [docs/abaqus_user_path.md](docs/abaqus_user_path.md) — the Abaqus user's route

Contracts and objects:
- [docs/program_definition.md](docs/program_definition.md) — what this program is
- [docs/problem_setting.md](docs/problem_setting.md) — the problem it solves
- [docs/input_objects.md](docs/input_objects.md) / [docs/output_objects.md](docs/output_objects.md)
- [docs/blackbox_order2_contract.md](docs/blackbox_order2_contract.md) — Path B, order ≥ 2
- [docs/residual_provider_contract.md](docs/residual_provider_contract.md) — Path C
- [docs/oti_execution_algorithm.md](docs/oti_execution_algorithm.md) — the order-by-order OTI loop
- [docs/simple_config_contract.md](docs/simple_config_contract.md) — the `resasm.yml` contract
- [docs/verification_contract.md](docs/verification_contract.md) — what is checked, and what a check does not prove
- [docs/privacy_contract.md](docs/privacy_contract.md) — private vs public
- [docs/which_path_should_i_use.md](docs/which_path_should_i_use.md)
- [docs/glossary.md](docs/glossary.md)

Also: [QUICKSTART_USER.md](QUICKSTART_USER.md) (one page),
[residual_core/README.md](residual_core/README.md), [STATUS.md](STATUS.md).

## Status

Honest accounting. Full detail in [STATUS.md](STATUS.md).

### What is verified

**Assembly of R is verified offline to machine precision** for the C3D8 solid
backends: divergence-theorem patch tests, a linear-stress patch test, frame
objectivity (5e-16), and an FD-checked tangent. `truss2`, `beam2`,
`nonlinear_spring1` and `nonlinear_bar1` assemble and are covered by the green
offline suite. The requirements engine, registry auto-selection, parser and
neutral-model IO are all exercised by tests.

**Sensitivities are verified** for the scalar proof backends (`nonlinear_spring1`,
`nonlinear_bar1`) and for the Python/black-box `resasm.yml` paths. STATUS.md records
the OTILib suite passing in WSL (11 passed, `RUN_OTILIB_TESTS=1`, so a skip would
have been a failure), including the order-2 recovery-factor proof.

### The gap you need to know about

**You cannot yet OTI-differentiate a C3D8+UMAT model through Path A.**

| backend | assemble R | OTI-differentiate R |
|---|---|---|
| `solid_c3d8_finite_strain`, `solid_c3d8_small_strain` | yes (verified) | **no** |
| `stress_driven_c3d8` | yes (verified) | **no** — σ is a frozen exported field |
| `truss2`, `beam2` | yes | **no** |
| `nonlinear_spring1`, `nonlinear_bar1` | yes | **yes** (proven in WSL) |

**Why:** the element and material kernels allocate numpy **float** arrays. For
example `residual_core/core/voigt.py::isotropic_D` builds
`np.zeros((6, 6), dtype=float)`, so a hypercomplex number cannot be stored in it and
the call fails. This is a concrete, fixable engineering gap — **not** a physics
limit.

The tool states this itself rather than failing late: `resasm init-assembly` /
`inspect-model` / `check` print a **Capability** block
(`assemble R: yes` / `OTI-differentiate R: NO`, naming the blocking backend), gated
in `resasm_user/recipe.py::sensitivity_capability`. `stress_driven_c3d8` is doubly
blocked — its σ is a frozen exported field, so it carries no parameter dependence to
differentiate.

**So, plainly, for an Abaqus user today:** you can **assemble and verify R** from
exported ingredients — that is the verified part. For *sensitivities* the routes are
Path B (black-box), or an **OTI-transformed UMAT** (the companion UMAT
source-transformation project) plugged into this assembler. That integration is the
intended design and the **next step — it is not a shipped feature.**

### Other real limits

- **Elements:** C3D8 is the only supported Abaqus solid. C3D8R / C3D20R / C3D4 are
  planned. `shell_placeholder` is a contract only, not runnable. `uel_direct` is a
  skeleton.
- **Loads:** only `*Cload` is assembled. `*Dsload` is **parsed but not applied**
  (`residual_core/core/loads.py`). Body forces likewise.
- **Constraints:** only Dirichlet / symmetry BCs are applied. `*Equation` / MPCs are
  **parsed but not applied** (`residual_core/core/constraints.py`).
- **The real CP UMAT** does not compile under gfortran (Cray pointers + an ifort
  `trace()` kind mismatch). It needs **Intel ifort + Abaqus**.
- **The comparison of our assembled R against Abaqus reaction forces is built but
  has never been executed** — Abaqus is not installed in this environment. ODB
  export, real-UMAT replay and the RF comparison are ready to run, not run. They are
  marked pending, never passing.
- **The finite-strain consistent tangent** is analytic but approximate, pending the
  Abaqus `AMATRX` comparison.
- **On Windows the Python residual path cannot run at all**: `oti_global.solve_python`
  hard-requires OTILib, which does not build natively on Windows. Use WSL.

The framework is formulation-agnostic *by architecture*. Each formulation becomes
supported when a backend satisfying the contract is registered and verified — not
before.

## Licensing

The redistributable core builds only against `sources/permissive/` (MIT / BSD-3).
`sources/copyleft/` and `sources/license-unknown/` are reference-only and are never
copied or linked. OTILib is GPLv3 and stays external. See `sources/LICENSING.md`.
