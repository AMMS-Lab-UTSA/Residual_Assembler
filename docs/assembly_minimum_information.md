# Assembly — the minimum information

What each term of the residual costs you, in inputs. Nothing here is a support
claim beyond what the code does.

```
    R(u, a)  =  F_internal(u, a, q)  -  F_external(a, t)  +  F_constraints(u, t)
    R_e      =  ∫_Ωe  Bᵀ σ(u, a, q) dΩ   -   f_e^ext
```

## The table

| Residual term | Required ingredient |
|---|---|
| internal force | mesh + element formulation + material + solution |
| external force | loads / stimuli |
| constraints | BCs / prescribed DOFs |
| state-dependent stress update | state/history + time increment |
| tangent | material tangent, **or** assembled tangent, **or** solver-exported tangent |
| **OTI RHS** | **an OTI-safe material/formulation, or a black-box coefficient provider** |

The last row is where the current gap lives. Read it twice.

---

## Row by row

### 1. Internal force — mesh + formulation + material + solution

`F_internal = Σ_e ∫ Bᵀσ dΩ`. Assembled by
`core/assembler.py` calling `Formulation.eval_element` per element and scattering
the result. The assembler contains no element mathematics and no material
knowledge.

| ingredient | concrete key | if missing |
|---|---|---|
| mesh | `model.inp` / `model.json` (positional arg to every `resasm` command) | nothing assembles; `requirements` reports `mesh` |
| element formulation | inferred from element type; override with `Config.formulation_policy = {C3D8: solid_c3d8_finite_strain}` | the element is **skipped** (`diag['skipped_no_formulation']`); `resasm inspect` says "backend missing" |
| material | `Config.subroutine` / `--subroutine`, or `Config.material_backend`, or a `section` dict | `material-replay` reports `material_model` missing. **Or**: skip the material entirely and use stress-driven mode. |
| solution `u` | `ResidualProblem.set_solution(U)` (API only — there is no CLI flag) | defaults to **zeros**. In stress-driven mode via `ResidualProblem` this is harmless (`options['config'] = 'small'`, so `u` never enters `F_internal`); in material-replay it means you assembled the residual at `u = 0`. |

**σ can come from two places** — this is the key structural fact:

| source of σ | mode | backend | needs a material? |
|---|---|---|---|
| an exported field (Abaqus' own IP stress) | `stress-driven` | `stress_driven_c3d8` | **no** |
| a material update from `u` | `material-replay` | `solid_c3d8_small_strain` / `solid_c3d8_finite_strain` | yes |

### 2. External force — loads / stimuli

`F_external`, assembled by `core/loads.py::external_force`.

| ingredient | concrete key | if missing |
|---|---|---|
| concentrated nodal loads | `*Cload` in the `.inp`; `cloads` in a neutral model | `F_external = 0` — **silently**. The residual then contains internal force only. |

**Not assembled, ever, today:**

| stimulus | what happens |
|---|---|
| `*Dsload` (distributed / pressure) | parsed into `AbaqusModel.dsloads`, **never applied** |
| body forces | not parsed, not applied |
| `*Amplitude` | not a parsed keyword (recorded in `unsupported_keywords`). `Boundary.amplitude` / `Cload.amplitude` hold the *name*, but `Assembler.assemble` calls `external_force(model, dm, load_factor)` with **no amplitude map**, so amplitude scaling never happens. |
| contact | not parsed, not applied |

If your model has any of these, `F_external` is incomplete and `R` will not
reproduce the solver's residual. The tool will not warn you.

### 3. Constraints — BCs / prescribed DOFs

| ingredient | concrete key | if missing |
|---|---|---|
| Dirichlet BCs | `*Boundary` (numeric DOF ranges, and `XSYMM`/`YSYMM`/`ZSYMM`/`XASYMM`/`YASYMM`/`ZASYMM`/`ENCASTRE`/`PINNED`) | every DOF is free; `||R_free||` then wrongly includes reaction DOFs, and no reaction is reported |

Constraints enter as a **partition**, not as an assembled force
(`core/constraints.py::partition` → `free_mask`, `prescribed_idx`). The residual
at the prescribed DOFs *is* the reaction (`Assembler.split`).

**Not applied:** linear `*Equation` / MPC constraints. They are parsed onto
`model.equations` and then ignored. Prescribed *values* are resolved
(`dirichlet_dofs` → `{dof: value}`) but nothing injects them into `u` — you must
supply a `u` that already satisfies them.

### 4. State-dependent stress update — state/history + time increment

Only needed when the material is history-dependent (plasticity, crystal
plasticity, viscoplasticity).

| ingredient | concrete key | if missing |
|---|---|---|
| previous state | per-element per-IP `STATEV` (count from `*Depvar`); managed by `core/state_manager.py` | `material-replay` reports `state_prev` missing. If you force it, the state slab is **zeros** — correct for an elastic law, wrong for a history-dependent one. |
| time increment | `assemble(time=(t_step, t_total), dtime=Δt)` | defaults to `time=(0,0)`, `dtime=0`. A rate-dependent law then sees a zero increment. `requirements` lists `time_increments`, but `ui/wizard.py::_availability` always reports it as present — **it does not verify that a real increment sequence was supplied**. |

The discipline (`residual_core/docs/limitations.md` §4): a history-dependent
replay must march the **entire** increment sequence, propagating state. It may
never jump to the final step and "solve" for it.

### 5. Tangent — three possible sources

`T = ∂R/∂u`. `core/results.py::TangentSource`:

| source | how you get it | which backends |
|---|---|---|
| `BACKEND_ASSEMBLED` | `assemble(compute_tangent=True)` — sums element `k_e` | `solid_c3d8_small_strain` (analytic, FD-verified to 1.6e-16), `solid_c3d8_finite_strain` (analytic, *approximate* — the `DDSDDE → AMATRX` mapping is open), `truss2`, `beam2`, `nonlinear_spring1`, `nonlinear_bar1` |
| solver-exported | `resasm.yml`: `tangent: {type: file, file: tangent.npz}` | any — the framework only needs the matrix |
| `UNAVAILABLE` | — | `stress_driven_c3d8` returns `k_e = None`. A stress-driven assembly has **no tangent**, by construction. |

If no element contributed a tangent, `ResidualProblem.result()` does not
fabricate one: it returns `TangentSource.UNAVAILABLE` with the reason attached
(`diag['tangent_contributions'] == 0`).

### 6. OTI RHS — the gap

This is the row that does not currently close for a real Abaqus solid model.

To produce `R^(p)` (the order-`p` hypercomplex residual coefficients) through
Path A, the residual must be **evaluated with OTI numbers**:
`core/oti_rhs_provider.py` seeds every parameter simultaneously, calls the same
`Formulation.eval_element`, and reads the coefficients off the returned residual.
That only works if the element kernel is written in **generic arithmetic**.

`Formulation.oti_differentiable` (in `formulations/base.py`) declares this. It is
`False` by default and `True` for exactly two backends.

## Verified capability matrix

| backend | assemble `R`? | OTI-differentiate `R`? |
|---|---|---|
| `solid_c3d8_small_strain` / `solid_c3d8_finite_strain` | ✅ verified offline to machine precision (divergence-theorem + linear-stress patch tests) | ❌ **NO** |
| `stress_driven_c3d8` (`R` from exported IP stress + `U`) | ✅ verified | ❌ NO (σ is a fixed exported field; `∂R/∂a = 0` by construction) |
| `truss2`, `beam2` | ✅ verified (`EA/L`, `PL³/3EI`) | ❌ NO |
| `nonlinear_spring1`, `nonlinear_bar1` | ✅ | ✅ **YES** — proven in WSL (`tests/framework/test_otilib_fe_sensitivity.py`) |
| `uel_direct` | adapter exists, needs a UEL routine | ❌ NO |
| `shell_placeholder` | contract only, NOT runnable | — |

### Why the solids are not OTI-differentiable

Not physics. **Float arrays.**

The element and material kernels allocate numpy `float` arrays and write scalars
into them:

- `core/voigt.py::isotropic_D` builds `D = np.zeros((6,6), dtype=float)`;
- `formulations/solid_c3d8_finite_strain.py` builds `sigma_ip = np.zeros((8,6))`
  and `D_ip = np.zeros((8,6,6))`;
- every solid `eval_element` starts with `u_e = np.asarray(dofs, dtype=float)`;
- `core/assembler.py` scatters with `np.add.at(R, edofs, np.asarray(r_e, float))`.

A float array cannot hold a hypercomplex number. Handed one, it either **raises**
(`ValueError: setting an array element with a sequence` — the failure mode
documented on `Formulation.oti_differentiable`) or **silently truncates** it to
its real part. Both destroy the imaginary directions, which are the derivative.

So, today:

- you **can** assemble `R` for a C3D8 model;
- you **cannot** yet OTI-overload the parameters through it.

Making the solid kernels generic-arithmetic / OTI-safe is a concrete, tractable
engineering task. It has not been done.

### The two routes that do work today

| route | how | what you get |
|---|---|---|
| **Path B — black-box** | your solver returns the order-`p` residual **coefficients** in a response file; the framework only solves `T U^(p) = -R^(p)` | full sensitivities, model stays private. Contract: [blackbox_order2_contract.md](blackbox_order2_contract.md). **Return coefficients, not derivatives** — at order ≥ 2 the difference is a factor of `κ!`. |
| **OTI-transformed UMAT** | source-transform the UMAT so it computes in OTI arithmetic (the companion UMAT source-transformation project), then assemble through this framework | **OTI-transformed UMAT + this assembler = an OTI-differentiable assembled residual.** This is the intended integration. It is the **next step, not a shipped feature** — no such UMAT is wired in this tree. |

Note that the OTI-transformed-UMAT route also requires making the *formulation*
kernel OTI-safe (the material returns a hypercomplex σ, which then has to survive
`sigma_ip`, `Bᵀσ`, and the assembler's scatter). Both halves are the same fix.

---

## Inferable vs. what you must supply

| item | inferable? | source of the inference | consequence if you rely on the default |
|---|---|---|---|
| DOF numbering | ✅ | `DofManager.for_model` — per-node union of the formulations' `dof_types` | none; mixed truss/beam/solid models work |
| element → formulation | ✅ | element type → registry (`default_formulation_policy`, `_policy_for_mode`) | unmapped element types (`C3D8R`, `C3D4`, …) are skipped, not assembled |
| element → material | ✅ | `*Solid Section, elset=…, material=…` | — |
| built-in material binding | ✅ | non-user `*Material` → `isotropic_elastic` | falls back to `[E, nu] = [200000, 0.3]` if the `.inp` has no constants |
| UMAT material | ❌ | — | `material-replay` blocked; use stress-driven |
| parameter *names* | ✅ | `section` keys and `PROPS[i]` (`_default_parameters`) | if the model has none, it falls back to `["parameter_1"]` → zero RHS |
| parameters to seed | ❌ | your choice | — |
| loads / BCs | ✅ (parsed) | `*Cload`, `*Boundary` | only `*Cload` is *assembled*; `*Dsload`, amplitudes, MPC are not |
| solution `u` | ⚠ zeros | — | material-replay residual evaluated at `u = 0` |
| exported stress field | ❌ | a prior solver run | stress-driven blocked |
| state / `STATEV` | ⚠ zeros | — | wrong for any history-dependent material |
| time / `dtime` | ⚠ `(0,0)`, `0.0` | — | wrong for any rate-dependent material |
| tangent | ✅ when the backend has one | element `k_e` | `UNAVAILABLE` in stress-driven mode — not an error, a fact |
| OTI-safety | ❌ | — | `oti_differentiable = False` → no Path-A sensitivities |

---

## Minimum viable input

### To assemble and verify `R` (the verified path)

```bash
resasm assemble model.inp --mode stress-driven --fields fields.json
resasm verify   model.inp --fields fields.json
```

Needs: **mesh + an exported integration-point stress field.** That is all — no
material, no PROPS, no state, no tangent (`core/requirements.py`:
`stress-driven` = `mesh`, `dof_field`, `element_field`; the `dof_field` check is
satisfied by the zero default).

Field format (`scripts/extract_odb_fields.py` writes it):

```json
{ "stress_ip": { "1": [[s11,s22,s33,s12,s13,s23], "... 8 IP rows"] } }
```

### To assemble `R` from a material (replay)

```bash
resasm requirements model.inp --mode material-replay
resasm assemble    model.inp --mode material-replay --subroutine umat.for --tangent
```

Needs: `mesh`, `solution_history`, `material_model`, `material_parameters`,
`state_prev`, `time_increments`. For a `*User Material`, `ResidualProblem`
refuses to fake it:

> `material-replay for user material 'X' needs the compiled UMAT (Intel ifort +
> Abaqus). Use mode='stress-driven' with an exported field for an offline
> residual.`

### To get parameter sensitivities

```bash
resasm sensitivity model.json --param mat.k --order 2 --backend otilib
```

Needs an `oti_differentiable` formulation (`nonlinear_spring1`,
`nonlinear_bar1`) **and** OTILib (Linux/WSL only — `scripts/run_otilib_tests_wsl.sh`).
For anything else, use Path B (black-box) via `resasm.yml`:

```yaml
residual:
  type: executable
  command: ./my_solver --request {request} --response {response}
tangent:
  type: file
  file: tangent.npz
parameters: {C11: 168400.0, C12: 121400.0}
solution: {file: converged_u.npy}
sensitivity: {order: 2, backend: otilib}
```

---

## See also

- [residual_assembly_recipe.md](residual_assembly_recipe.md) — the full ingredient list.
- [abaqus_user_path.md](abaqus_user_path.md) — the Abaqus route, and what it cannot reconstruct.
- [blackbox_order2_contract.md](blackbox_order2_contract.md) — coefficients, not derivatives.
- `residual_core/docs/limitations.md` · `STATUS.md`
