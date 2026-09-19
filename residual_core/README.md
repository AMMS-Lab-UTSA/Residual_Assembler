# residual_core: the formulation-agnostic residual assembly framework

This page describes the `residual_core` package for developers and reviewers:
its two interfaces, its layout, the registered backends and how to run the
offline verification. Users who want to compute sensitivities should start at
[START_HERE.md](../START_HERE.md) and [README.md](../README.md).

`residual_core` assembles the **global finite-element residual**

```
R(y, q, a, t) = F_internal(y, q, a, t) - F_external(t) + F_constraints(y, t)
```

outside the FE solver, for arbitrary element formulations and materials. The
**core assembler is formulation-agnostic**: it makes no assumption about
displacement-only DOFs, element topology, mechanics, stress and strain, UMATs
or crystal plasticity. It only collects element residuals and scatters them.
Everything physics-specific lives behind two interfaces, a **Formulation**
(element weak form) and a **Material** (constitutive update), and is
discovered through a backend **registry**.

> **The framework is formulation-agnostic by architecture. Each formulation
> becomes supported when a backend satisfying the formulation contract is
> registered and verified.**

The same package holds the engines behind the sensitivity commands:
`replay/` replays a compiled UMAT-OTI provider against an Abaqus analysis and
solves for total-history parameter sensitivities (`resasm request`,
`resasm history`, `resasm replay`); `algebra/` connects OTILib for the
hypercomplex sensitivity paths.

**Crystal plasticity is only one example backend, not the organising
principle.** Verified backends today:

- **stress-driven C3D8 solid**: residual from an externally exported stress field;
- **C3D8 material replay**: small- and finite-strain solids with isotropic
  elasticity, the bounded compressible neo-Hookean law, or a UMAT (crystal
  plasticity is one example, kept intact);
- **truss/bar** (`truss2`) and **3D beam/frame** (`beam2`): proofs that the core
  is not C3D8- or CP-specific, exercising heterogeneous DOF sets;
- **nonlinear spring and bar** (`nonlinear_spring1`, `nonlinear_bar1`): the
  proof backends for OTI differentiation through the assembler;
- **UEL-direct** residual adapter (skeleton);
- **shell** backend *contract* (placeholder: declared and documented, not runnable).

What is supported and how it is verified: [STATUS.md](../STATUS.md). The
backend-level limitations: [docs/limitations.md](docs/limitations.md).

## Quick start

```python
from residual_core import ResidualProblem

p = ResidualProblem.from_abaqus("model.inp")
p.inspect()                              # elements, materials, reachable modes, missing data
R = p.assemble(mode="stress-driven")     # after attaching an exported field
```

```bash
resasm inspect  model.inp                                # what can be assembled, what is missing
resasm assemble model.inp --mode stress-driven --fields fields.json
resasm doctor   model.inp                                # per-mode readiness + minimum missing input
# (use `python -m residual_core.ui.cli ...` if no console script is installed)
```

See [docs/user_interface.md](docs/user_interface.md) and
[docs/minimal_input_contract.md](docs/minimal_input_contract.md). Every
command is described in [docs/CLI_GUIDE.md](../docs/CLI_GUIDE.md).

## The two interfaces

```python
# formulations/base.py  — one element family's weak form
eval_element(element_id, element_type, coords, dofs, solution_state,
             material_state, properties, time, dtime, fields, options)
    -> (element_residual, element_tangent, updated_state, diagnostics)

# materials/base.py  — one constitutive law (the UMAT contract, generalized)
evaluate(kinematics, state_prev, binding, time, dtime, fields, options)
    -> (stress, tangent, state_new, diagnostics)
```

`core/assembler.py::Assembler` loops over elements, dispatches each to its
bound `Formulation`, and scatters; it contains **no** element or material
mathematics. The evidence that it is genuinely agnostic:
[tests/framework/test_assembler.py](../tests/framework/test_assembler.py)
drives the CP example through the generic assembler and reproduces the
verified C3D8 kernel **bitwise** in both offline modes.

## Layout

```
residual_core/
  core/           formulation-agnostic engine
    assembler.py      loop elements -> eval_element -> scatter -> loads -> reactions
    model.py          neutral Model + from_abaqus() converter
    dof_manager.py    heterogeneous global DOF numbering (per-node DOF sets)
    state_manager.py  per-IP material state across increments (history-safe)
    constraints.py    Dirichlet BCs / free-prescribed split (equations parsed)
    loads.py          F_external (concentrated loads)
    registry.py       generic BackendSpec + Registry
    requirements.py   per-mode minimum-input engine (data minimization)
    diagnostics.py    model inspector (auto backend selection + reachable modes)
    verification.py   generic Levels 0-7 (zero-field, rigid-body, FD-tangent, reactions)
    where_it_went_wrong.py  nine-layer failure attribution (residual_core.diagnose)
  formulations/   element weak forms (Formulation backends) + registry.py
    base.py, c3d8_kernel.py (proven numerics),
    solid_c3d8_small_strain.py, solid_c3d8_finite_strain.py (home of the CP backend),
    stress_driven_adapter.py (Mode 1), uel_adapter.py (Mode 3),
    truss2.py, beam2.py, nonlinear_spring1.py, nonlinear_bar1.py,
    shell_base.py + shell_placeholder.py
  materials/      constitutive backends (Material) + registry.py
    base.py, elastic_adapter.py, neo_hookean.py, umat_adapter.py,
    crystal_plasticity_adapter.py (one example backend, behind the interface)
  replay/         compiled-provider replay and total-history sensitivities
                  (bounded J2 engine, history engine, ODB export to .npz)
  algebra/        OTILib adapter (production) and Dual1 (legacy first-order smoke test)
  io/             abaqus_inp_parser.py, abaqus_odb_export.py, neutral_model_io.py
  ui/             ResidualProblem facade, resasm CLI, config, runnable examples (no physics)
  app/            optional Streamlit front end (a thin layer over the CLI)
  examples/       manifest.py (classifies every .inp by residual backend)
  docs/           architecture, contracts, user_interface, minimal_input_contract,
                  adding_a_formulation, adding_a_material, verification_strategy, limitations
tests/
  framework/      the generic core is formulation-agnostic (Abaqus-free):
                  assembler, truss2, beam2, mixed-model dispatch
  cp_c3d8_umat/   the C3D8+UMAT CP backend's verification (one backend)
  replay_history/ the history replay engine
```

## The assembly modes

| Mode | Backend(s) | Field source | Runs offline? |
|------|-----------|--------------|---------------|
| **1 stress-driven** | `stress_driven_c3d8` | element field from an Abaqus or other export | yes (verified) |
| **2 material-replay** | `solid_c3d8_*` + a `materials/*` backend | material `evaluate()` at each IP | yes with `isotropic_elastic` and `compressible_neo_hookean`; the Oxford CP UMAT needs ifort + Abaqus |
| **3 UEL / direct** | `uel_direct` | a UEL-like routine's RHS/AMATRX | yes (mock); sign convention `R = -RHS` documented and verified |
| **4 formulation** | `truss2`, `beam2`, `nonlinear_spring1`, `nonlinear_bar1` | geometry + section properties | yes (verified: `EA/L`, `PL³/3EI`) |

Registered backends: formulations `solid_c3d8_finite_strain`,
`solid_c3d8_small_strain`, `stress_driven_c3d8`, `truss2`, `beam2`,
`nonlinear_spring1`, `nonlinear_bar1`, `shell_placeholder` (contract only) and
`uel_direct`; materials `isotropic_elastic`, `compressible_neo_hookean`, `umat`
and `crystal_plasticity`. The model inspector auto-selects a backend per
element type; the requirements engine reports the *minimum* missing input for
a mode. Adding backends: [docs/adding_a_formulation.md](docs/adding_a_formulation.md),
[docs/adding_a_material.md](docs/adding_a_material.md).

## Run the offline verification (no Abaqus needed)

```bash
python tests/framework/test_assembler.py                       # generic core == verified kernel (both modes)
python tests/framework/test_truss2_backend.py                  # bar backend: axial force EA/L + FD tangent
python tests/framework/test_beam2_backend.py                   # beam backend: cantilever PL^3/3EI + FD tangent
python tests/framework/test_mixed_model_dispatch.py            # physics-blind dispatch: truss + beam + solid in one model
python tests/framework/test_dof_manager_mixed.py               # heterogeneous DOF edge cases (union per node)
python tests/framework/test_requirements_negative.py           # missing-data reports (minimum missing input only)
python tests/framework/test_public_api.py                      # ResidualProblem facade only (no backend imports)
python tests/framework/test_neutral_io.py                      # solver-neutral JSON round-trip + resasm inspect model.json
python -m residual_core.ui.examples                            # runnable, self-checking API examples
python -m residual_core.examples.generate_minimal             # (re)generate the minimal CLI example assets
python residual_core/formulations/c3d8_kernel.py               # patch tests + FD tangent (machine precision)
python tests/cp_c3d8_umat/tangent_fd_check/run_checks.py       # CP backend, tabulated
python residual_core/stress_driven_residual.py --selftest      # global equilibrium (uniform + linear)
python -m residual_core.materials.umat_adapter                 # UMAT<->formulation contract (mock UMAT)
python residual_core/umat_adapter_fortran/umat_replay.py --dry-run   # UMAT driver plumbing
python residual_core/examples/manifest.py                      # regenerate the backend manifest
```

All of these run on Python 3 and NumPy. The full offline suite is
`python -m pytest -q -m "not abaqus and not arc and not network"` from the
repository root.

## Example inventory (by residual backend)

[example_manifest.md](example_manifest.md) and `.json` classify all **105**
`.inp` files: **A** standard solid + UMAT (6), **B** standard solid + built-in
(1), **C** UEL/direct (29), **D** coupled-field (0), **E** unsupported element
(13), **F** reference-only by licence (56). The four priority C3D8 CP jobs are
group A, backend `solid_c3d8_finite_strain`, status *supported*.

## Verification and status

[STATUS.md](../STATUS.md) records what is verified and how, including the
framework's original success criteria and the crystal-plasticity backend's
criteria. [docs/verification_strategy.md](docs/verification_strategy.md)
describes the Level 0 to 7 ladder and which levels each backend exercises.
Comparisons against Abaqus (tested: 2021.HF5) run where Abaqus is licensed:
the replay engines check replayed stress, state and reactions against the ODB
at every integration point and increment. Everything Abaqus-independent is
verified offline.

## Licensing

The redistributable core builds only against `sources/permissive/` (MIT and
BSD-3). `sources/copyleft/` (AGPL) and `sources/license-unknown/` are
reference-only: all 56 files are manifest group **F**, never copied or linked.
See [sources/LICENSING.md](../sources/LICENSING.md).
