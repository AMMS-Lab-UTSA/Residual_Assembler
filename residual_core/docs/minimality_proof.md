# Minimality Proof

Why the framework is **model-agnostic by architecture** and asks for the
**minimum** data — argued from the code, not asserted.

## 1. The core has no physics

`core/assembler.py::Assembler.assemble` contains no element or material
mathematics. Its per-element body is: resolve the bound formulation key, gather
the element's DOFs/coords/state/props, call `form.eval_element(...)`, and scatter
the returned residual/tangent. It never imports a concrete formulation or
material. Proof it works unchanged across unrelated physics:

- stress-driven C3D8, elastic C3D8, crystal-plasticity C3D8, **truss**, **beam**,
  and a **mixed** truss+beam+solid model all run through the same assembler
  (`tests/framework/test_assembler.py`, `test_truss2_backend.py`,
  `test_beam2_backend.py`, `test_mixed_model_dispatch.py`).

If the core embedded any physics, one of these would require special-casing. None
does.

## 2. Heterogeneous DOFs, not displacement-only

`core/dof_manager.py::DofManager.for_model` assigns each node the **union** of the
`dof_types` of the formulations touching it — 3 for a truss node, 6 for a beam
node, coexisting in one dense global numbering
(`tests/framework/test_dof_manager_mixed.py`). No step assumes all nodes share the
same DOFs.

## 3. Backends *declare*; the framework *reasons*

Every backend advertises a `BackendSpec` (`core/registry.py`): element types, DOF
types, supported modes, per-mode required/optional inputs, material-interface
need, state requirements, tangent support, verification status. The inspector
(`core/diagnostics.py`) and requirements engine (`core/requirements.py`) reason
purely over these declarations — never over physics. Adding a backend needs no
core change (`docs/adding_a_formulation.md`, `docs/adding_a_material.md`).

## 4. Minimum data, per mode

`core/requirements.py` reports the **single** minimum missing input for a chosen
mode, not a generic checklist (`tests/framework/test_requirements_negative.py`):

| Mode | Minimum inputs |
|---|---|
| 1 stress-driven | mesh + DOF field + exported element field |
| 2 material-replay | mesh + material + PROPS + history + state + dtime |
| 3 direct/UEL | mesh + DOF layout + callable UEL |
| 4 formulation | backend + section properties |

The material model is **not** required for Mode 1: the residual is assembled from
the exported field alone. That is the sharpest statement of minimality — the
assembly needs the *field that enters the weak form*, nothing more.

## 5. Crystal plasticity is one backend

CP lives entirely in `formulations/solid_c3d8_finite_strain.py` +
`materials/crystal_plasticity_adapter.py`, registered alongside six other
formulation backends and two other materials. Remove it and the core, the truss,
the beam, the stress-driven path, and the inspector all still work. It is a row in
`docs/support_matrix.md`, not the table.

## 6. External generality (the zoo)

The verification zoo (`tests/verification_zoo/`) maps external UMAT/UEL examples
across viscoelasticity, J2 plasticity, damage, hyperelasticity, cohesive, and
coupled chemo-mechanics. Every one reduces, for assembly purposes, to the **same**
stress-driven or direct-residual path already exercised internally — confirming
the agnostic claim beyond crystal plasticity, honestly (most are `planned` or
`reference-only`, each with a stated next step).

## 7. The sensitivity backend is pluggable too

Producing $\mathbf{R}^{(p)}$ is a **declared capability**, not core physics. A
`SensitivityRHSProvider` (`core/rhs_provider.py`) is selected by name; the
production one is `OtiLibRHSProvider` (OTILib, arbitrary order, all parameters
seeded simultaneously, solved order-by-order with lower-order $\mathbf{U}$ injected
before advancing). `DualNumberRHSProvider` (Dual1) remains only as a legacy
first-order smoke test. Remove either and the assembly core is unchanged; a missing
OTILib is reported cleanly rather than hidden behind Dual1.

## Conclusion

Agnosticism is structural: no physics in the core, heterogeneous DOFs, declared
capabilities, and per-mode minimum inputs. The claim the framework makes is the
honest one — *formulation-agnostic by architecture; a formulation becomes
supported when a backend satisfying the contract is registered and verified.*
