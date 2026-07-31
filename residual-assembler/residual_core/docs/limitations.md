# Limitations and Scope

This document states, plainly and honestly, what the framework does **not** do,
what it needs to support a new formulation, and the concrete gaps in the current
tree. It exists so the framework's guarantees are not overread. Facts only.

---

## 1. A UMAT alone is NOT enough to reconstruct a residual

An Abaqus UMAT returns, at a single material (integration) point, the stress, the
material tangent `DDSDDE`, and the updated state `STATEV`. It does **not** return
the element weak form: the shape functions, the B-operator, the integration rule,
or how those point quantities are integrated into element nodal forces. Abaqus
performs that element assembly internally and hides it. The residual
`f_int,e = integral_v B^T sigma dv` therefore cannot be recovered from a UMAT
alone — the element formulation (`formulations/`) is a separate, required piece.

## 2. A material model gives stress/tangent/state update, not element assembly

Restated at the level of the `Material` contract (`materials/base.py`):
`evaluate(...)` answers the constitutive question only — `(stress, tangent,
state_new)` at a point, in a declared stress/tangent measure. It knows nothing
about elements, shape functions, or assembly; that is the `Formulation`'s job
(`formulations/base.py`). A material is necessary but not sufficient for a
residual.

## 3. To support a NEW formulation the framework needs one of four things

Supporting a new element family means providing the weak form in one of these
forms (each maps to a backend mode in `docs/architecture.md`):

1. the **element weak form** (shape functions + integration + B-operator + how
   the field is integrated), implemented as a `Formulation` — e.g. the
   `solid_c3d8_*` backends;
2. a **known standard element implementation** whose kernel can be reused (the
   verified `c3d8_kernel` is the example);
3. **exported integration-point fields** (stress-driven, Mode 1) — the field
   entering the weak form is supplied externally
   (`formulations/stress_driven_adapter.py`); or
4. a **UEL-like direct-residual routine** returning element RHS/AMATRX directly
   (Mode 3, `formulations/uel_adapter.py`).

Without one of these, the framework cannot assemble that formulation's residual.

## 4. History-dependent materials require the full loading history

Plasticity (including crystal plasticity) is history-dependent: the stress at an
increment depends on the accumulated state, not just the current strain. A replay
or residual evaluation must **march the entire loading path increment by
increment, propagating state** — it may never jump to the final step and "solve"
for it. This discipline is enforced by `core/state_manager.py::StateManager`,
which separates committed (increment-start) state from trial state and only
commits an accepted increment (`get`/`set_trial`/`commit`/`rollback`). The
standalone UMAT replay (`residual_core/umat_adapter_fortran/umat_replay.py`)
likewise persists `STATEV` + the `/UMPS/` common block across increments and warns
when output frames are too sparse for an exact history replay.

## 5. OTI/HYPAD sensitivity comes LATER; the large model is NOT started

- **OTI/HYPAD** (operator-overloading / hypercomplex automatic sensitivity) is
  **not implemented**. It will build on this residual framework later, but nothing
  in the current tree provides it.
- The **150-parameter model is NOT started**.
- The framework does **not** claim arbitrary-formulation support. It is
  formulation-agnostic *by architecture*: a formulation is supported once a
  backend satisfying the contract is registered **and verified**. The currently
  registered backends (`formulations.default_formulations()` /
  `build_formulation_registry()`) are: `solid_c3d8_finite_strain`,
  `solid_c3d8_small_strain`, `stress_driven_c3d8`, `truss2`, `beam2`,
  `shell_placeholder` (contract only, not runnable), and (optionally)
  `uel_direct`. Any other formulation is unsupported until built per item 3.

---

## Current concrete gaps (as of this tree)

- **Continuum solids: only C3D8 has a verified backend.** The default policy maps
  `C3D8 -> solid_c3d8_finite_strain`; C3D8R / C3D20R / C3D4 (and U1/cohesive) have
  **no default continuum backend** — they are kept in the model but not assembled
  by the default solid path ("planned", not implemented). Non-solid families now
  have proof backends: `truss2` (bar, `T3D2`/`T2D2`) and `beam2` (3D frame,
  `B31`/`B33`), which are small-strain linear elements — proofs of agnosticism,
  not production elements. A `shell_placeholder` declares the shell contract but
  is not runnable.
- **CP material runs offline only with a mock UMAT.** The Mode-2 finite-strain
  backend needs a `Material` with `kinematic_input='deformation_gradient'`.
  `materials/` ships `elastic_adapter.py::IsotropicElastic` (runnable),
  `umat_adapter.py::UmatAdapter`, and `crystal_plasticity_adapter.py::CrystalPlasticityAdapter`;
  `default_materials()` registers `['crystal_plasticity','isotropic_elastic','umat']`.
  `CrystalPlasticityAdapter` IS wired into the `Assembler` as a `Material` (through
  `SolidC3D8FiniteStrain`). Its `backend='python'` path runs offline with an
  injected mock UMAT (isotropic elasticity — NOT crystal plasticity); its
  `backend='fortran'` path shells out to the standalone driver
  (`residual_core/umat_adapter_fortran/`) and raises a clear ifort+Abaqus error if
  the real binary is not built (it refuses the mock binary and never fabricates a
  stress). So Levels 6/7 DO have an in-tree material to drive — they only need the
  ifort-built `umat_driver` to produce genuine crystal-plasticity stress.
- **The real CP UMAT needs ifort + Abaqus.** The Grilli `umat.for` does not compile
  under gfortran (Cray-pointer/`target` twin arrays + an ifort `trace()` kind
  mismatch); it requires Intel **ifort + Abaqus (MKL)**. Until then only the
  **mock** elastic UMAT runs, so the material-replay comparison is against the
  mock, not crystal plasticity.
- **Abaqus C3D8 IP ordering is unvalidated offline.** `ABAQUS_C3D8_GAUSS` uses the
  standard lexicographic (xi1-fastest) order; the internal `sigma_ip[k] <->
  points[k]` pairing is verified, but the match to Abaqus' export index is
  provably not testable with the uniform/self-sampled stress fields used offline.
  It must be confirmed against the first real ODB with a spatially-varying stress
  before trusting Mode-1 on a non-uniform state.
- **Heterogeneous mechanical DOFs work; coupled physics DOFs are not implemented.**
  `core/dof_manager.py` builds per-node DOF sets from the formulations touching
  each node, so translational (`UX,UY,UZ`) and rotational (`RX,RY,RZ`) DOFs coexist
  in one model (`beam2`, and the mixed-dispatch test). However, no *non-mechanical*
  DOF (temperature, pressure, concentration) has a formulation or constitutive
  wiring yet — coupled/thermal problems are designed-for by the DOF machinery but
  not implemented.
- **Constraints beyond Dirichlet are parsed-not-applied.**
  `core/constraints.py` resolves `*Boundary` (numeric ranges + symmetry keywords
  XSYMM/YSYMM/ZSYMM/ENCASTRE/PINNED/...) into the free/prescribed partition. Linear
  `*Equation` / MPC constraints are carried on the model (`model.equations`) but
  **not applied** in the residual. Distributed loads (`*Dsload`) and body forces
  are likewise parsed elsewhere but not assembled — `core/loads.py` handles only
  concentrated `*Cload`.

---

## Documentation-vs-code note

The refactored code is the source of truth. Some in-repo prose predates the
refactor and describes the earlier flat layout: `residual_core/README.md`,
`residual_core/CONTRACT.md`, and `STATUS.md` still refer to a top-level
`c3d8_residual.py` (now `formulations/c3d8_kernel.py`),
`extract_abaqus_fields.py` (now `io/abaqus_odb_export.py`), and test paths
`tests/c3d8_tangent_fd_check/` etc. (now under `tests/cp_c3d8_umat/`), and they
state a deliberately narrow "C3D8 + UMAT only, no UEL path, no generic elastic
MVP" scope that the refactor has since broadened (a formulation-agnostic core, a
UEL adapter, and a runnable elastic material). Where those documents and the code
disagree, the code and these `docs/` files govern.
