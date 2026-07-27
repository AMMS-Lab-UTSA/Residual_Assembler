# The Residual Assembly Recipe

## The reframe

The old framing asked the user to *provide the residual*: "write
`residual(u, params)`". For a real solver like Abaqus that is **not possible** —
the solver owns the global residual internally and never exposes it. There is no
call, no output request, no ODB field that hands you `R`.

The correct framing:

> **The user should not provide `R`. The user should provide enough ingredients
> for us to build `R`.**
>
> Abaqus hides `R`, but `R` is not magic. It is assembled from element residuals,
> materials, solution fields, loads and constraints. So we ask for *those*, not
> for `R`.

The residual we assemble:

```
    R(u, a)  =  F_internal(u, a, q)  -  F_external(a, t)  +  F_constraints(u, t)

    R_e      =  ∫_Ωe  Bᵀ σ(u, a, q) dΩ   -   f_e^ext
```

`u` = solution field · `a` = parameters · `q` = state/history · `t` = time.

The central object is therefore a **residual assembly recipe**: the list of
ingredients from which `R` can be rebuilt outside the solver.

### What the code actually assembles today

`residual_core/core/assembler.py::Assembler.assemble` computes

```
    R = Σ_e scatter( formulation.eval_element(...) )  -  F_external
```

and `Assembler.split(R)` partitions `R` into free DOFs and prescribed DOFs (the
reaction). Constraints are handled **by partition, not by an assembled force
term**: `core/constraints.py` resolves `*Boundary` (numeric + symmetry keywords)
into the free/prescribed split. `F_constraints` above is the general form; in the
current code there is no assembled constraint-force contribution — linear
`*Equation`/MPC constraints are carried on the model and **not applied**
(`core/constraints.py` module docstring; `residual_core/docs/limitations.md`).

`F_external` is concentrated nodal load only (`core/loads.py`: `*Cload`).

---

## Analogy with the UMAT source-transformation project

That project's contract was a transformation JSON. This project's contract is an
assembly recipe. Same idea, different ingredients:

```
Old UMAT project:        source + seed + output + target + promote + replace
New Residual_Assembler:  mesh + formulation + material + fields + stimuli + parameters
```

Both say: *give me the pieces and the rules; I will produce the differentiable
artifact.*

---

## Three paths, in priority order

| Path | What the user supplies | Where it lives in the code | Status |
|---|---|---|---|
| **A — assembly from ingredients** (primary) | mesh, element formulation, material, solution field, state, loads/BCs, parameters, tangent source | `residual_core` — `resasm inspect / requirements / assemble`, `ResidualProblem` API | Assembles `R`. **Not** OTI-differentiable for solids (see below). |
| **B — black-box** | a private executable that returns residual *coefficients* per order | `resasm_user` — `resasm.yml` with `residual.type: executable`; contract in [blackbox_order2_contract.md](blackbox_order2_contract.md) | Working. The framework never sees the model. |
| **C — direct residual** (shortcut/toy) | a Python `residual(u, params, state, time)` returning the global `R` | `resasm.yml` with `residual.type: python` | Working. **Not** the main story — a real solver cannot do this. |

> **Naming warning.** The older docs in this folder
> ([which_path_should_i_use.md](which_path_should_i_use.md),
> [residual_provider_contract.md](residual_provider_contract.md)) letter these
> differently: there, "Path A" is the *Python global residual* and "Path C" is the
> *black-box executable*. This page uses the ingredient-first lettering above.
> When in doubt, refer to the `residual.type` value (`python` / `executable` /
> `element`) or to the assembly mode name — those are unambiguous and appear in
> the code.

### Which surface implements Path A

Path A is implemented in **`residual_core`**, driven by the `resasm inspect /
requirements / assemble / verify` commands and the `ResidualProblem` Python API.

It is **not** reachable from `resasm.yml`. `residual.type: element` is accepted by
`resasm_user/config.py`, but `resasm_user/runner.py` and `resasm_user/checks.py`
treat `element` **identically to `python`**: they load `residual.module` and call
a global `residual(u, params, state, time)`. No element-wise assembly happens on
that path. (The prose in `which_path_should_i_use.md` and `user_input_contract.md`
that says the framework "assembles element-by-element" for `residual.type: element`
overstates what the runner does.)

---

## The recipe objects

Eleven objects. For each: what it is, where it comes from, whether it can be
inferred, and what happens if it is missing.

> **There is no `recipe.yml` parser.** The recipe is a *specification* — the
> checklist of ingredients. Today each object is expressed through one of four
> concrete surfaces:
>
> | surface | file / call |
> |---|---|
> | model file | `model.inp` (Abaqus) or `model.json` (`residual_core/io/neutral_model_io.py`, schema `resasm-neutral-model/1`) |
> | advanced config | `residual_core/ui/config.py::Config` (`--config x.yml`): `mode`, `odb`, `subroutine`, `formulation_policy`, `material_backend`, `material_parameters`, `options` |
> | CLI flags | `resasm assemble MODEL --mode M --fields F --subroutine S --tangent --out R.npy` |
> | Python API | `ResidualProblem.from_abaqus / from_neutral / attach_results / attach_subroutine / set_solution / assemble / result / sensitivity_package` |

### 1. `mesh`

Nodes, elements, element types, connectivity, node/element sets.

| | |
|---|---|
| Source | `*Node` / `*Element` / `*Nset` / `*Elset` in the `.inp` (`io/abaqus_inp_parser.py`), or `nodes` / `elements` in a neutral `model.json` |
| Inferable | No. This is the one thing that must always be given. |
| Consumed by | `core/model.py::Model`, `core/dof_manager.py` |
| If missing | Nothing can be assembled. The requirements engine reports `mesh` missing for every mode (`core/requirements.py`). |

Limits: `*Include` is spliced; a non-identity `*Instance` transform raises
`NotImplementedError` (the parser refuses to silently drop it).

### 2. `formulation`

The element weak form: DOFs, shape functions, integration rule, B-operator, how
the field entering the weak form becomes nodal forces. Contract:
`formulations/base.py::Formulation.eval_element`.

| | |
|---|---|
| Source | a registered backend (registry key), chosen per element |
| Inferable | **Yes**, from the element type. `core/model.py::default_formulation_policy` maps `C3D8 → solid_c3d8_finite_strain`; `ui/wizard.py::_policy_for_mode` re-selects per requested mode (e.g. `C3D8 → stress_driven_c3d8` in stress-driven mode). Override with `Config.formulation_policy = {ELEMENT_TYPE: backend}`. |
| If missing | The element is kept in the model but **not assembled** (`assembler.py` counts it in `diag['skipped_no_formulation']`). `resasm inspect` reports "backend missing" and names the next step. |

Registered formulations (`formulations/registry.py`):

| backend | element types | DOFs/node | modes | status |
|---|---|---|---|---|
| `solid_c3d8_small_strain` | `C3D8` | UX,UY,UZ | material-replay | verified |
| `solid_c3d8_finite_strain` | `C3D8` | UX,UY,UZ | material-replay | verified |
| `stress_driven_c3d8` | `C3D8` | UX,UY,UZ | stress-driven | verified |
| `truss2` | `T3D2`,`T2D2`,… | UX,UY,UZ | formulation | implemented/simple |
| `beam2` | `B31`,`B33`,… | UX,UY,UZ,RX,RY,RZ | formulation | implemented/simple |
| `nonlinear_spring1` | `SPRING1`,… | U | formulation | verified |
| `nonlinear_bar1` | `NLBAR2`,… | U | formulation | verified |
| `uel_direct` | `U1` (configurable) | user-declared | direct-residual | skeleton (needs a UEL routine) |
| `shell_placeholder` | `S3`,`S4`,`S4R`,… | 6 | **none** | contract-only — `eval_element` raises |

Not implemented: `C3D8R`, `C3D20R`, `C3D4`, `C3D10` — planned, no backend. They
parse, they inspect, they do not assemble.

### 3. `material`

The constitutive law at a point: `(stress, tangent, state_new)`. Contract:
`materials/base.py::Material.evaluate`. A material is **necessary but not
sufficient** — it does not know shape functions or assembly
(`residual_core/docs/limitations.md` §1–2).

| | |
|---|---|
| Source | `*Material` / `*User Material` / `*Depvar` in the `.inp`; or `MaterialBinding` in the API; or `section` dict in a neutral model |
| Inferable | Partly. A non-user material is bound to `isotropic_elastic` automatically (`ui/wizard.py::_bind_materials_for_replay`, default `[E, nu] = [200000, 0.3]` if no constants). A `*User Material` is **not** inferable: it needs the UMAT. |
| Registered | `isotropic_elastic` (runnable), `umat`, `crystal_plasticity` (`materials/registry.py`) |
| If missing | `material-replay` mode reports `material_model` missing. In stress-driven mode the material is **not needed at all** — the stress field replaces it. |

`kinematic_input` must match the formulation: `small_strain` (`strain`,`dstrain`)
for `solid_c3d8_small_strain`; `deformation_gradient` (`F0`,`F1`) for
`solid_c3d8_finite_strain`.

### 4. `fields` — solution field and exported element fields

Two distinct things.

**Solution field `u`** (nodal DOFs):

| | |
|---|---|
| Source | `ResidualProblem.set_solution(U)`, or `assemble(U=...)`, or solved internally by `solve_newton()` |
| Inferable | Defaults to **zeros** if not given (`ui/wizard.py::_run`). There is no CLI flag for `u`. |
| Note | In stress-driven mode through `ResidualProblem`, `options['config']` is hard-set to `'small'` (reference configuration), so `u` **does not enter** `F_internal` at all. The current-configuration (nlgeom) assembly of an exported stress *does* use `u` — see `residual_core/stress_driven_residual.py --mode finite`. |

**Exported element field** (the field entering the weak form — for solids, the
integration-point Cauchy stress):

| | |
|---|---|
| Source | `--fields fields.json` / `attach_results(...)` |
| Shape | `{"stress_ip": {"<eid>": [[s11,s22,s33,s12,s13,s23] × n_ip]}}`, Abaqus Voigt order |
| Produced by | `scripts/extract_odb_fields.py` (run under Abaqus Python) |
| Inferable | No. It comes from a prior solver run. |
| If missing | `stress-driven` mode reports `element_field` missing and refuses to assemble. |

Caveat: when `--fields` points at a JSON **path**, `attach_results` keeps only
`stress_ip`; exported `reactions` / `displacements` / `statev` in the same file
are **dropped**. Pass a dict through the API (or use
`scripts/compare_residuals.py`, which reads `reactions` itself) to use them.

### 5. `stimuli` — loads and boundary conditions

| | |
|---|---|
| Source | `*Cload`, `*Boundary` in the `.inp`; `cloads` / `boundaries` in a neutral model |
| Inferable | No — but they are parsed automatically from the model file. |
| Assembled | **Only concentrated `*Cload`** (`core/loads.py::external_force`). |
| **Not assembled** | `*Dsload` (parsed into `AbaqusModel.dsloads`, never applied), body forces, contact, amplitudes (see below). |
| If missing | `F_external = 0`. The residual is then internal force only — physically wrong for a loaded model, and silently so. |

Amplitudes: `*Amplitude` is not a parsed keyword (it lands in
`unsupported_keywords`). `Boundary.amplitude` / `Cload.amplitude` record the
*name* only, and `Assembler.assemble` calls `loads.external_force(model, dm,
load_factor)` without an amplitude map — so amplitude scaling is **never
applied**. Use `load_factor` yourself if you need it.

### 6. `state` — history variables

| | |
|---|---|
| Source | UMAT `STATEV` / ODB `SDV`; `*Depvar` gives the count |
| Managed by | `core/state_manager.py::StateManager` (committed vs trial state; `get`/`set_trial`/`commit`/`rollback`) |
| Inferable | Only trivially: zeros. That is correct for a stateless material and **wrong** for a history-dependent one. |
| If missing | `material-replay` reports `state_prev` missing. History-dependent materials must march the full increment sequence — never jump to the final step (`residual_core/docs/limitations.md` §4). |
| Not needed for | stress-driven mode (the stress is already the answer). |

### 7. `time` — increment sequence

| | |
|---|---|
| Source | the caller: `assemble(time=(step_t, total_t), dtime=Δt)` |
| Inferable | Defaults to `time=(0.0, 0.0)`, `dtime=0.0`. |
| If missing | Rate-dependent / viscoplastic materials get `dtime = 0` and are wrong. The requirements engine lists `time_increments` for `material-replay` but `ui/wizard.py::_availability` reports it as always available — it does **not** check that a real increment sequence was supplied. |

### 8. `constraints`

| | |
|---|---|
| Applied | Dirichlet `*Boundary`: numeric DOF ranges and symmetry keywords `XSYMM/YSYMM/ZSYMM/XASYMM/YASYMM/ZASYMM/ENCASTRE/PINNED` (`core/constraints.py::_SYMM_DOF`) |
| Effect | free/prescribed partition; the assembled residual at prescribed DOFs is the reaction |
| **Not applied** | linear `*Equation` / MPC (parsed into `model.equations`, never used); contact; tie constraints |
| If missing | Everything is free. `||R_free||` then includes what should have been reaction DOFs, and no reaction can be reported. |

Prescribed **values** are resolved (`dirichlet_dofs` returns `{dof: value}`) but
nothing injects them into `u` — you must supply a `u` that already satisfies them.

### 9. `parameters` — the design variables `a`

| | |
|---|---|
| Source | `MaterialBinding.section` dict (`{"E": 210000, "nu": 0.3}`), or `constants` (PROPS), or `Config.material_parameters = {matname: [PROPS…]}` |
| Naming | `"<material>.<key>"`, e.g. `spring.k`; `ResidualProblem._default_parameters()` enumerates section keys and `PROPS[i]` entries |
| Used by | `resasm sensitivity MODEL --param mat.key` / `sensitivity_package(parameters=[...])` |
| Inferable | The *names* are, from the material's section/constants. Which ones you want to differentiate is not. |
| If missing | `_default_parameters()` falls back to `["parameter_1"]`, which matches nothing and yields a zero RHS. Name them explicitly. |

Backends that declare `sensitivity_parameters`: `truss2` (`E`,`A`),
`nonlinear_spring1` (`k`,`f`), `nonlinear_bar1` (`k`).

### 10. `tangent` — `T = ∂R/∂u`

Three sources, in `core/results.py::TangentSource`:

| source | how |
|---|---|
| backend-assembled | `assemble(compute_tangent=True)` — sums element `k_e`. Available for `solid_c3d8_*` (analytic), `truss2`, `beam2`, `nonlinear_*`. |
| solver-exported | supply `T` yourself (`resasm.yml`: `tangent.type: file`, `tangent.file: tangent.npz`) |
| unavailable | `stress_driven_c3d8` returns `k_e = None` — a stress-driven assembly has **no** tangent. `ResidualProblem.result()` marks it `TangentSource.UNAVAILABLE` with the reason. |

If no element contributed a tangent, `result()` does not fabricate one — it says
so (`diag['tangent_contributions'] == 0`).

Open item: mapping a finite-strain UMAT `DDSDDE` to Abaqus' `AMATRX` (objective
rate + geometric split) is **not closed** — the finite-strain element tangent is
approximate pending the Abaqus comparison (declared in the backend's own
`tangent_support` string, and in `STATUS.md` criterion 7).

### 11. `sensitivity` — the request

| | |
|---|---|
| Fields | parameters, derivative order `q`, hypercomplex backend |
| Backends | `otilib` (production, arbitrary order — `core/oti_rhs_provider.py`) or `dual1` (legacy, order 1 only — `core/rhs_provider.py`) |
| Method | seed all parameters simultaneously, then for `p = 1..q`: evaluate `R*` with OTI scalars, extract `R^(p)`, solve `T U^(p) = -R^(p)`, inject `U^(p)` back into `u*` |
| Requires | an **OTI-safe** formulation, i.e. `Formulation.oti_differentiable == True` |
| Platform | OTILib does not build on Windows. Verified in WSL (`scripts/run_otilib_tests_wsl.sh`); on Windows the OTI tests skip cleanly. |

**`oti_differentiable` is `True` for exactly two backends today:**
`nonlinear_spring1` and `nonlinear_bar1` (grep `oti_differentiable` — the base
class default is `False`).

Everything else — including every C3D8 solid backend — assembles `R` but cannot
carry an OTI number through its kernel, because those kernels use numpy **float**
arrays (`np.zeros((6,6), dtype=float)` in `core/voigt.py::isotropic_D`,
`np.zeros((8,6))` for the IP stress slab, `np.asarray(dofs, dtype=float)` in the
element entry points, `np.asarray(r_e, float)` in the assembler's scatter). A
hypercomplex scalar handed to those arrays is either rejected outright or
truncated to its real part — either way the imaginary directions are destroyed.
This is a **fixable engineering gap (make the kernels generic-arithmetic), not a
physics limitation**; it is what `formulations/base.py::oti_differentiable`
documents.

---

## The full recipe

Every object, written out. This is the *specification* view — see the surface
table above for where each line actually lives today.

```yaml
# ResidualAssemblyRecipe — the ingredients of R(u, a)
mesh:
  file: model.inp                 # or model.json (resasm-neutral-model/1)
  # nodes, elements, element types, connectivity, sets

formulation:
  policy:                         # -> Config.formulation_policy
    C3D8: solid_c3d8_finite_strain
  # inferable from element type; override only to force a different backend

material:
  CPuranium:
    backend: umat                 # -> Config.material_backend
    source: umat.for              # -> Config.subroutine / --subroutine
    parameters: [ ... ]           # PROPS -> Config.material_parameters
    n_state_vars: 125             # from *Depvar

fields:
  solution: U.npy                 # -> ResidualProblem.set_solution(U)  (API only)
  exported: fields.json           # -> --fields ; {"stress_ip": {eid: [[6] x n_ip]}}

stimuli:
  loads:      from_model          # *Cload only  (assembled)
  # *Dsload / body force / contact: parsed or ignored, NEVER assembled
  boundaries: from_model          # *Boundary (Dirichlet + symmetry) -> partition

state:
  previous: statev.json           # STATEV / SDV per element per IP
  history:  full_increment_sequence   # required for history-dependent materials

time:
  step:   [t_step, t_total]       # -> assemble(time=...)
  dtime:  1.0e-3                  # -> assemble(dtime=...)

constraints:
  dirichlet: applied              # from *Boundary
  equations: NOT APPLIED          # *Equation / MPC parsed, never assembled

parameters:                       # the design variables a
  - CPuranium.PROPS[3]
  - CPuranium.PROPS[7]

tangent:
  source: backend_assembled       # | solver_exported (tangent.npz) | unavailable

sensitivity:
  order: 2
  backend: otilib                 # requires an oti_differentiable formulation
```

## The minimal recipe

Almost all of the above is inferred. The stress-driven verification path — the
one that is actually verified today — needs exactly this:

```yaml
mesh:
  file: model.inp
fields:
  exported: fields.json           # integration-point stress
```

which is one command:

```bash
resasm assemble model.inp --mode stress-driven --fields fields.json
```

Everything else is inferred: the formulation (`C3D8 → stress_driven_c3d8`), the
DOF numbering, the loads and BCs (parsed from the `.inp`), the solution field
(defaults to zeros — and does not enter the reference-configuration assembly),
the material (not needed), the state (not needed), the tangent (none — a
stress-driven assembly has no tangent).

## Inferable vs must-be-supplied

| object | inferable? | from what |
|---|---|---|
| `mesh` | ❌ never | — |
| `formulation` | ✅ | element type → registry (`default_formulation_policy`, `_policy_for_mode`) |
| DOF layout | ✅ | union of the formulations' `dof_types` per node (`DofManager.for_model`) |
| `material` (built-in) | ✅ | non-user `*Material` → `isotropic_elastic` |
| `material` (UMAT) | ❌ | needs the subroutine, or use stress-driven |
| `parameters` (names) | ✅ | material `section` keys / `PROPS[i]` |
| `parameters` (which to seed) | ❌ | you must choose |
| `stimuli` | ✅ (parsed) | `*Cload` / `*Boundary` in the model file — but only `*Cload` is assembled |
| `fields.solution` | ⚠ defaults to zeros | supply it, or `solve_newton()` |
| `fields.exported` | ❌ | a prior solver run |
| `state` | ⚠ defaults to zeros | wrong for history-dependent materials |
| `time` | ⚠ defaults to `(0,0)`, `dtime=0` | wrong for rate-dependent materials |
| `tangent` | ✅ when the backend has one | otherwise export it or do without |
| `sensitivity` | ❌ | order + parameters + backend are a request, not data |

## Progressive disclosure — the tool tells you what is missing

The requirements engine (`core/requirements.py`) reports the **single minimum
missing item**, never a generic checklist:

```bash
resasm inspect      model.inp            # what can be assembled, by which backend
resasm requirements model.inp --mode stress-driven
resasm doctor       model.inp            # per-mode readiness + a config template
resasm backends                          # every backend's declared capabilities
```

Minimum inputs per mode, verbatim from `MODE_REQUIREMENTS`:

| mode | minimum inputs |
|---|---|
| `stress-driven` | `mesh`, `dof_field`, `element_field` |
| `material-replay` | `mesh`, `solution_history`, `material_model`, `material_parameters`, `state_prev`, `time_increments` |
| `direct-residual` | `mesh`, `element_dofs`, `uel_routine` |
| `formulation` | `formulation_backend`, `section_properties` |

## See also

- [assembly_minimum_information.md](assembly_minimum_information.md) — what each
  residual term costs you in inputs, and where the OTI gap is.
- [abaqus_user_path.md](abaqus_user_path.md) — the Abaqus-specific route and its
  exact limits.
- [blackbox_order2_contract.md](blackbox_order2_contract.md) — Path B.
- `residual_core/docs/limitations.md` — the backend-level limitation list.
- `STATUS.md` — what is verified vs what is built-but-not-executed.
