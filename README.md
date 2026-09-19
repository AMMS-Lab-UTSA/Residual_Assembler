# Residual_Assembler

[![CI](https://github.com/AMMS-Lab-UTSA/Residual_Assembler/actions/workflows/ci.yml/badge.svg)](https://github.com/AMMS-Lab-UTSA/Residual_Assembler/actions/workflows/ci.yml)
[![License: GPL-3.0-only](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

Residual_Assembler computes **parameter sensitivities of a converged
finite-element analysis without re-running it**. From the saved analysis
(`Analysis.inp` + `Analysis.odb`), a compiled material provider
(`OTI_UMAT.obj` + `Mapping.json`, built by the companion
[UMAT-OTI](https://github.com/AMMS-Lab-UTSA/UMAT_source_transformation)) and a
`sensitivity_request.json`, it replays the material at every integration point
and increment, assembles the residual `R`, its tangent `K` and `dR/dp`, solves

    K du/dp = -dR/dp

and writes the requested sensitivities (`sensitivity_results.json`,
`sensitivity_tables.csv`, `run_report.txt`; full fields under `private/`). The
material's source code never reaches the person running the analysis.

```bash
resasm request --model Analysis.inp --odb Analysis.odb \
    --material OTI_UMAT.obj --request sensitivity_request.json --out results
```

## Quick start

```bash
git clone https://github.com/AMMS-Lab-UTSA/Residual_Assembler.git
git clone https://github.com/AMMS-Lab-UTSA/UMAT_source_transformation.git
python3 -m venv .venv && . .venv/bin/activate
pip install -e "./Residual_Assembler[gui,yaml,test]" "./UMAT_source_transformation[test]"
resasm --help
```

Linux with Python 3.10+ and `gfortran`. Abaqus (tested: 2021.HF5) is needed
only to read an ODB; the sensitivity computation itself does not call it.

**New here?** Install with [docs/INSTALL.md](docs/INSTALL.md), then work
through the seven [worked examples](examples/README.md), from a one-line
residual to full-size Abaqus models. The [command-line guide](docs/CLI_GUIDE.md)
and the [GUI guide](docs/GUI_GUIDE.md) cover every entry point with real
output; [docs/USAGE_REPORT.md](docs/USAGE_REPORT.md) summarises the program as
it works now.

## Worked examples

Seven complete examples, each with a walkthrough: what it shows, the exact
commands, the equivalent GUI steps, the measured output and how it is checked
independently ([examples/README.md](examples/README.md)).

| # | Example | What it shows | Abaqus? |
| --- | --- | --- | --- |
| 1 | [The smallest sensitivity](examples/user_config_minimal/WALKTHROUGH.md) | A residual you write yourself (a cubic spring); first and second derivatives against the closed form | No |
| 2 | [Stress-driven C3D8 assembly](residual_core/examples/minimal_c3d8_stress_driven/WALKTHROUGH.md) | The residual of one element from a supplied stress field, against the analytic face tractions | No |
| 3 | [The four-file request on one element](examples/presentation_request/WALKTHROUGH.md) | `Analysis.inp` + `Analysis.odb` + `OTI_UMAT.obj` + `sensitivity_request.json` to sensitivities; J2 plasticity against the uniaxial closed form | To read the ODB |
| 4 | [History replay of a J2 beam, offline](examples/replay_history/WALKTHROUGH.md) | Ten plastic increments: ODB parity, Abaqus finite differences, the homogeneity identity, weighted shares, re-equilibration | No (export committed) |
| 5 | [Full-size cantilevers](examples/cantilevers/WALKTHROUGH.md) | J2 (1,536 C3D8, 40 increments) and FCC crystal plasticity (384 C3D8, 25 increments, 10 parameters); full-field sensitivities | Once, to run the analyses |
| 6 | [Provider-to-sensitivity pipeline](examples/bounded_j2_c3d8/WALKTHROUGH.md) | Builds the compiled J2 material, solves a cyclic one-element history, verifies every derivative against whole-model finite differences | No |
| 7 | [Finite-strain neo-Hookean C3D8](examples/finite_strain_c3d8/WALKTHROUGH.md) | Finite-strain assembly with the exact tangent and parameter sensitivities, against nonlinear re-solves | No |

Measured on the two cantilevers: the replayed stress, state and reactions match
the ODB at every integration point; the J2 model runs in about 10 s (25 s when
every increment is first re-equilibrated to double precision). The
sensitivities agree with whole-model central finite differences of the original
UMAT to 6.6e-8 (E), 5.8e-8 (nu), 7.9e-9 (initial yield stress) and 5.4e-6 (H),
and to at most 6.7e-7 for all ten FCC parameters. A check that needs no finite
differences holds as well: both models are homogeneous of degree one in their
stress-dimensioned parameters, so at every increment the parameter-weighted
sensitivities of each stress and reaction sum to the value itself, and those of
each displacement sum to zero. The re-equilibrated results satisfy this to
1e-12 (J2) and 1.2e-13 (FCC). Method and full results:
[docs/REPLAY_HISTORY.md](docs/REPLAY_HISTORY.md),
[docs/VERIFICATION_RECORD.md](docs/VERIFICATION_RECORD.md).

## How it works

**You do not have to provide R. You provide the ingredients, and
Residual_Assembler builds R.** Abaqus keeps its global residual internal, so
"write your residual function" is not something an Abaqus user can do. But R is
assembled from element residuals, materials, solution fields, loads and
constraints, and those are what the program asks for:

```
R(u, a) = F_internal(u, a, q) - F_external(a, t) + F_constraints(u, t)

R_e     = ∫_Ωe B^T σ(u, a, q) dΩ - f_e^ext
```

Supported, each with its own verified workflow: C3D8 stress-driven assembly;
total-history sensitivities of small-strain C3D8 analyses with any UMAT-OTI
provider (`resasm history`, which `resasm request` uses automatically for every
model outside its bounded single-material example); and finite-strain
neo-Hookean sensitivities. Finite-strain plasticity, other element types,
several steps or materials and distributed loads are refused with a named
reason rather than approximated.

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
  formulation backend: yes
  solution field (U / U+rotation / T): yes
  stress / resultant field: no
  one-step deck in scope: yes

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

### The GUI

There is an optional Streamlit front end over the same commands:

```
pip install -e ".[gui]"
streamlit run scripts/app.py
```

The first tab, **Sensitivity Request**, is the collaborator screen: point at `OTI_UMAT.obj` (its `Mapping.json` beside it), the
saved `Analysis.inp` and `Analysis.odb`, tick the parameters, choose the output
and the region, and press **Solve** — it runs `resasm request` and shows the
full-field result ([docs/GUI.md](docs/GUI.md), screenshots in
`docs/screenshots/`). The other tabs — Start here, Model, Requirements,
Assemble, Sensitivity, Job, Backends, Advanced Replay — mirror the CLI one-for-one. The
GUI is deliberately thin: every button builds an argv list and calls
`residual_core.ui.cli.main` in-process, then shows the command it ran, the real
exit code and the captured output verbatim, so it cannot show a number the CLI
would refuse to produce (`tests/framework/test_gui_is_a_thin_cli_front_end.py`
pins that). The interactive `resasm init` wizard stays CLI-only.

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

Requires Python ≥ 3.9 (3.10 with UMAT-OTI), NumPy and SciPy. `PyYAML` is
optional (`pip install -e .[yaml]`); a minimal YAML reader is bundled. The
Streamlit GUI is also optional (`pip install -e ".[gui]"`); nothing in
`residual_core` imports it. For the collaborator workflow install UMAT-OTI
beside it; [docs/USAGE_REPORT.md](docs/USAGE_REPORT.md) gives the full recipe
and the clean-install gate that checks it.

Path A (assembly, inspection, verification) needs nothing else from PyPI, but
the offline test suite reads one **external source submodule**. A fresh clone
does not contain it, so run:

```
./scripts/init_permissive_sources.sh    # fetches the permissive test dependency
```

That fetches only `sources/permissive/` (MIT / BSD-3). Copyleft and
license-unknown submodules are marked `update = none` and are never fetched by
setup. Without the bootstrap the tests that need it skip with a message naming
the missing file and this command. The offline suite is
`pytest -q -m "not abaqus and not arc and not network"`, with OTILib and the
UMAT-OTI checkout beside this one; its result from clean clones of both
repositories is in [docs/evidence/final_clean_clone.md](docs/evidence/final_clean_clone.md). See
[sources/SUBMODULES.md](sources/SUBMODULES.md) for the tier policy, the pinned
commits, and the clean-clone verification procedure.

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

Start here:
- [docs/INSTALL.md](docs/INSTALL.md) — installation from an empty machine, the checks that prove it works, troubleshooting
- [examples/README.md](examples/README.md) — seven worked examples
- [docs/CLI_GUIDE.md](docs/CLI_GUIDE.md) — every command, its options, real output and exit codes
- [docs/GUI_GUIDE.md](docs/GUI_GUIDE.md) — a step-by-step walkthrough of every screen
- [docs/USAGE_REPORT.md](docs/USAGE_REPORT.md) — the program as it works now, in one place
- [docs/REQUEST_INTERFACE.md](docs/REQUEST_INTERFACE.md) — the four-input `resasm request` interface
- [docs/REPLAY_HISTORY.md](docs/REPLAY_HISTORY.md) — `resasm history`: the mathematics, tolerances and supported deck subset
- [docs/GUI.md](docs/GUI.md) — the GUI's design, with screenshots
- [docs/VERIFICATION_RECORD.md](docs/VERIFICATION_RECORD.md) — every verified result, its command, reference and measured value
- [docs/CONNECTED_WORKFLOW.md](docs/CONNECTED_WORKFLOW.md) — the provider-to-sensitivity pipeline in one command
- [docs/COMPATIBILITY.md](docs/COMPATIBILITY.md) — the pinned UMAT-OTI version and the shared contract

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

### Sensitivities of C3D8 models with a UMAT

**The shipped route is the compiled OTI provider.** `resasm request` (bounded
single-material J2 engine) and `resasm history` (any provider, prescribed
displacements, long histories, sparse assembly) replay an OTI-transformed UMAT
built by UMAT-OTI and solve for full-field sensitivities; `resasm request`
hands a model to the history engine when it is outside the bounded scope and
says so. Verified against whole-model finite differences of the original UMAT
and against Abaqus reruns — see [docs/REPLAY_HISTORY.md](docs/REPLAY_HISTORY.md).

**Path A with the Python OTILib backend** (the `resasm.yml` recipe that
differentiates the built-in C3D8 kernels in pyoti) is still reported as
`OTI-differentiate R: NO` for the C3D8 backends by `init-assembly` /
`inspect-model` / `check`; that capability gate
(`resasm_user/recipe.py::sensitivity_capability`) is deliberate and stays until
that path is verified. `stress_driven_c3d8` carries a frozen exported stress
and therefore no parameter dependence to differentiate.

### Other real limits

- **Elements:** C3D8 is the only supported Abaqus solid. C3D8R / C3D20R / C3D4 are
  planned. `shell_placeholder` is a contract only, not runnable. `uel_direct` is a
  skeleton.
- **Loads:** only `*Cload` is assembled. `*Dsload` is **parsed but not applied**
  (`residual_core/core/loads.py`). Body forces likewise.
- **Constraints:** only Dirichlet / symmetry BCs are applied. `*Equation` / MPCs are
  **parsed but not applied** (`residual_core/core/constraints.py`).
- **Deck keywords that are not applied are named.** `resasm inspect` and
  `resasm doctor` list them under **Deck keywords present but NOT applied**
  (`*Dsload`, `*Dload`, `*Equation`, `*Amplitude`, `AMPLITUDE=` references and any
  keyword the reader leaves unread, such as `*Elastic`), and every command that
  reads a deck prints them as a `DeckKeywordNotApplied` warning. `resasm verify`
  refuses (exit 2) a deck whose loads or constraints are not all applied.
- **One step at a time on the general deck path.** `inspect`, `requirements`,
  `assemble` and `verify` apply every `*Cload` of a deck at once and every
  `*Boundary` as one list, and integrate stress-driven and small-strain elements
  over the reference configuration. A deck with several steps, or `NLGEOM=YES`
  outside a finite-strain backend, is reported as out of scope and refused, with
  the reason. `resasm history` replays a deck step by step.
- **Material replay uses the model's own constants.** A deck's `*Elastic` is not
  read, so material replay of such a deck is refused (`needs: material parameters
  (PROPS)`, with the reason) rather than run with invented values.
- **The real CP UMAT** does not compile under gfortran (Cray pointers + an ifort
  `trace()` kind mismatch). It needs **Intel ifort + Abaqus**.
- **Abaqus comparisons** run where Abaqus 2021 is licensed: the replay engines
  check replayed stress, state and reaction forces against the ODB at every
  integration point and increment (single-precision ODB tolerances are
  documented in each run report). Without Abaqus, an ODB cannot be exported and
  the command says so.
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
copied or linked. This framework is GPL-3.0-only, so OTILib's GPLv3 is
compatible; it stays an external dependency because it is a separate
product, not because of licensing. See `sources/LICENSING.md`.

### Equilibrium verification exit status

`resasm verify MODEL --fields FIELDS --atol 1e-6` checks the Euclidean norm of
all free residual entries against an absolute tolerance in the model's force
units. It exits 0 only for finite residuals within tolerance, 1 for failed
equilibrium, and 2 when verification cannot run (including a deck whose loads
or constraints the residual does not apply). Choose `--atol` for the units
and precision of the exported analysis. Assembled reactions are printed, but
are explicitly **not checked against a reference** by this command. The shipped
minimal stress-driven cube has stress without balancing loads and therefore
correctly fails equilibrium verification, while still demonstrating assembly.
