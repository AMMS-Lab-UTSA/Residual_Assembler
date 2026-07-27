# residual_core — Formulation-Agnostic Residual Assembly Framework

> Full methodology (inputs, outputs, the residual method, testing, and the OTIS UMAT
> connection): [../METHODOLOGY.md](../METHODOLOGY.md).

Assemble the **global finite-element residual**

```
R(y, q, a, t) = F_internal(y, q, a, t) - F_external(t) + F_constraints(y, t)
```

outside the FE solver, for arbitrary element formulations and materials. The
**core assembler is formulation-agnostic**: it makes no assumption about
displacement-only DOFs, element topology, mechanics, stress/strain, UMAT, or
crystal plasticity. It only knows how to collect element residuals and scatter
them. Everything physics-specific lives behind two interfaces — a **Formulation**
(element weak form) and a **Material** (constitutive update) — and is discovered
through a backend **registry**.

> **The framework is formulation-agnostic by architecture. Each formulation
> becomes supported when a backend satisfying the formulation contract is
> registered and verified.**

**Crystal plasticity is only one example backend, not the organizing principle.**
Verified example backends today:

- **stress-driven C3D8 solid** — residual from an externally exported stress field;
- **UMAT-driven C3D8 crystal plasticity** — one example backend (kept intact);
- **simple truss/bar** (`truss2`) and **simple 3D beam/frame** (`beam2`) — proofs
  the core is not C3D8/CP-specific, exercising heterogeneous DOF sets;
- **UEL-direct** residual adapter (skeleton);
- **shell** backend *contract* (placeholder — declared, documented, not yet runnable).

> Guardrails (unchanged): no OTI/HYPAD yet, no 150-parameter model, no claim of
> arbitrary-formulation support beyond registered+verified backends. See
> [docs/limitations.md](docs/limitations.md).

## Quick start

```python
from residual_core import ResidualProblem

p = ResidualProblem.from_abaqus("model.inp")
p.inspect()                              # elements, materials, reachable modes, missing data
R = p.assemble(mode="stress-driven")     # after attaching an exported field
```

```bash
resasm inspect  model.inp                            # what can be assembled, what's missing
resasm assemble model.inp --mode stress-driven --odb fields.json
resasm doctor   model.inp                            # per-mode readiness + minimum missing input
# (use `python -m residual_core.ui.cli ...` if no console script is installed)
```

See [docs/user_interface.md](docs/user_interface.md) and
[docs/minimal_input_contract.md](docs/minimal_input_contract.md).

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

The `core/assembler.py::Assembler` loops elements, dispatches each to its bound
`Formulation`, and scatters — containing **no** element or material mathematics.
Proof it is genuinely agnostic: [tests/framework/test_assembler.py](../tests/framework/test_assembler.py)
drives the CP example through the generic assembler and reproduces the verified
C3D8 kernel **bitwise** in both offline modes.

## Layout

```
residual_core/
  core/           formulation-AGNOSTIC engine
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
  formulations/   element weak forms (Formulation backends) + registry.py
    base.py, c3d8_kernel.py (proven numerics),
    solid_c3d8_small_strain.py, solid_c3d8_finite_strain.py (home of the CP backend),
    stress_driven_adapter.py (Mode 1), uel_adapter.py (Mode 3),
    truss2.py, beam2.py (proofs of agnosticism), shell_base.py + shell_placeholder.py
  materials/      constitutive backends (Material) + registry.py
    base.py, elastic_adapter.py (runnable), umat_adapter.py,
    crystal_plasticity_adapter.py (one example backend, behind the interface)
  io/             abaqus_inp_parser.py, abaqus_odb_export.py, neutral_model_io.py
  ui/             ResidualProblem facade, resasm CLI, config, runnable examples (no physics)
  examples/       manifest.py (classifies every .inp by residual backend)
  docs/           architecture, contracts, user_interface, minimal_input_contract,
                  adding_a_formulation, adding_a_material, verification_strategy, limitations
tests/
  framework/      the generic core is formulation-agnostic (Abaqus-free):
                  assembler, truss2, beam2, mixed-model dispatch
  cp_c3d8_umat/   the C3D8+UMAT CP backend's verification (ONE backend)
```

## The assembly modes

| Mode | Backend(s) | Field source | Runs offline? |
|------|-----------|--------------|---------------|
| **1 stress-driven** | `stress_driven_c3d8` | element field from an Abaqus/other export | yes (verified) |
| **2 material-replay** | `solid_c3d8_*` + a `materials/*` backend | material `evaluate()` at each IP | yes with `isotropic_elastic`; CP needs ifort+Abaqus |
| **3 UEL / direct** | `uel_direct` | a UEL-like routine's RHS/AMATRX | yes (mock); sign convention `R = -RHS` documented + verified |
| **4 formulation** | `truss2`, `beam2` | geometry + section properties | yes (verified: `EA/L`, `PL³/3EI`) |

Registered backends: formulations `solid_c3d8_finite_strain`,
`solid_c3d8_small_strain`, `stress_driven_c3d8`, `truss2`, `beam2`,
`shell_placeholder` (contract only), `uel_direct`; materials `isotropic_elastic`,
`umat`, `crystal_plasticity`. The model inspector auto-selects a backend per
element type; the requirements engine reports the *minimum* missing input for a
mode. New backends: [docs/adding_a_formulation.md](docs/adding_a_formulation.md),
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
All of the above pass today on Python 3 + numpy.

## Example inventory (by residual backend)

[example_manifest.md](example_manifest.md) / `.json` classify all **105** `.inp`:
**A** standard solid + UMAT (6), **B** standard solid + built-in (1), **C** UEL/direct (29),
**D** coupled-field (0), **E** unsupported element (13), **F** reference-only/license (56).
The 4 priority C3D8 CP jobs are group A, backend `solid_c3d8_finite_strain`, status *supported*.

## Verification & status

See [../STATUS.md](../STATUS.md) for the honest accounting (both the CP Step-10
criteria and the refactor's success criteria) and [docs/verification_strategy.md](docs/verification_strategy.md)
for the Level 0–7 ladder and which levels each backend exercises. Abaqus is not
installed here, so Level-5/6/7 (ODB comparison, real-UMAT replay) are delivered
ready-to-run and clearly marked pending; everything Abaqus-independent is verified.

## Licensing

Redistributable core builds only against `sources/permissive/` (MIT/BSD-3).
`sources/copyleft/` (AGPL) and `sources/license-unknown/` are reference-only —
all 56 are manifest group **F**, never copied or linked. See `sources/LICENSING.md`.
