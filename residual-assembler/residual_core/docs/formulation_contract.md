# Formulation Contract — the 11 items every backend must define

A *formulation* owns the finite-element weak form for one element family. The
formulation-agnostic core (`core/assembler.py`) will drive **any** backend that
subclasses `formulations/base.py::Formulation` and implements `eval_element`, but
a backend is not *trusted* until all 11 items below are pinned down. This
document states the contract and fills in what the C3D8 solid backends
(`formulations/solid_c3d8_small_strain.py`,
`formulations/solid_c3d8_finite_strain.py`, numerics in
`formulations/c3d8_kernel.py`) actually provide as the worked example.

Conventions referenced throughout (from `residual_core/CONTRACT.md`): Voigt order
`(11,22,33,12,13,23)` (Abaqus order); engineering shear
`eps_voigt = (e11,e22,e33, 2e12, 2e13, 2e23)`; element DOFs node-major.

---

## 1. Number and type of DOFs per node

Declared by the class attribute `dof_types` (and `dofs_per_node = len(dof_types)`).
The `DofManager` and the assembler use this to build/gather global DOF indices.

**C3D8 solid:** `dof_types = ("UX","UY","UZ")` -> 3 displacement DOFs/node, 24 per
element, node-major ordering
`[u1x,u1y,u1z, ..., u8x,u8y,u8z]`.

## 2. Element interpolation / shape functions

The trial/test interpolation used to map nodal DOFs to the field and its
gradient.

**C3D8 solid:** trilinear shape functions
`N_a(xi) = 1/8 (1+xi1*xi1a)(1+xi2*xi2a)(1+xi3*xi3a)`, provided by
`c3d8_kernel.shape_functions` / `shape_grad_natural`; node natural coordinates in
`ABAQUS_C3D8_NODES` (8x3), in Abaqus C3D8 connectivity order.

## 3. Integration rule

The quadrature points and weights, and (critically for pairing with exported
fields) their ordering.

**C3D8 solid:** full 2x2x2 Gauss, points at `+/-1/sqrt(3)`, all weights 1, exposed
as the named constant `c3d8_kernel.ABAQUS_C3D8_GAUSS` (`.points` (8,3),
`.weights` (8,)) in Abaqus integration-point (SDV/S) output order (xi1 fastest).
The *internal* pairing `sigma_ip[k] <-> points[k]` is verified; the *match to
Abaqus' export index* is the one convention not offline-verifiable (see item 11
and `docs/verification_strategy.md`).

## 4. Kinematic variables computed from the DOFs

What the formulation derives from `dofs` before calling the material / forming the
weak form.

**C3D8 solid:**
- small strain (`solid_c3d8_small_strain`): engineering-shear strain
  `eps = B0 @ u_e` and increment `deps = B0 @ (u_e - u_prev)`, using the
  reference-config B-matrix `b_matrix_reference`.
- finite strain (`solid_c3d8_finite_strain`): deformation gradients
  `F1 = I + du/dX` (current) and `F0 = I + du_prev/dX` (previous) per IP, via
  reference-config shape-function gradients (`_F_at`), passed to the material as
  `{"F0","F1","element","ip"}`.

## 5. State variables required at integration points

The per-IP history slab the formulation reads (previous) and returns (updated),
sized `(n_ip, n_state)`. `StateManager` holds committed vs trial copies.

**C3D8 solid:** `n_state` comes from the `MaterialBinding.n_state_vars` (Abaqus
NSTATV / `*Depvar`). The formulation allocates `state_new = (8, nstate)`, reads
`material_state[ip]` as `s_prev`, writes `s_new` back per IP. The elastic example
(`IsotropicElastic`) is stateless (`n_state_vars = 0`).

## 6. Material / update law interface

How the formulation obtains stress + tangent. Modes 2 backends call a `Material`
(`materials/base.py`); Mode 1 reads an exported field; Mode 3 gets RHS/AMATRX
directly.

**C3D8 solid (Mode 2):** requires a `MaterialBinding` in `properties`; calls
`material.evaluate(kin, s_prev, binding, time, dtime, fields, options)` and expects
`(stress_voigt(6,), D(6,6) or None, state_new, diag)`. Small-strain passes
`kinematic_input='small_strain'` data; finite-strain passes
`kinematic_input='deformation_gradient'` data. Raises a clear error if no material
binding is attached (directs the caller to `stress_driven_adapter` for Mode 1).

## 7. Residual contribution

The `+F_internal` element vector returned as the first tuple element.

**C3D8 solid:**
- small strain: `r = sum_k B0^T sigma_k detJ0_k w_k` (reference config).
- finite strain: `r = c3d8_kernel.element_internal_force_finite_strain(coords,
  u_e, sigma_ip)` = `sum_k B_spatial(x)^T sigma_cauchy_k detJ_current_k w_k`, i.e.
  Cauchy stress integrated over the current configuration `x = X + u` — the
  nlgeom=YES Abaqus-matching path.

## 8. Tangent contribution (if available)

The `(ndof,ndof)` element tangent, or `None`.

**C3D8 solid:**
- small strain: `K = sum_k B0^T D_k B0 detJ0_k w_k` (exact for the linear case).
- finite strain: `c3d8_kernel.element_tangent(coords, u_e, D_ip, sigma_ip,
  mode="finite")` = material term `B_spatial^T D B_spatial` + conventional
  geometric/initial-stress term `kron(G, I3)` with `G = dNdx sigma dNdx^T`.
  Returned only when `options['compute_tangent']` is set.
- **Open item:** mapping the finite-strain UMAT `DDSDDE` into the exact Abaqus
  `AMATRX` (objective-rate + geometric split) is deferred to the Abaqus
  comparison; the residual (Cauchy stress only) is exact, the geometric term and
  the small-strain material term are verified. `c3d8_kernel` also exposes
  `force_tangent_fixed_sigma` (the exact fixed-sigma force Jacobian), which is the
  correct analytic target of the finite-strain force FD check and is *distinct*
  from the geometric stiffness.

## 9. External load contribution (if supported)

Whether the backend contributes `F_external` (tractions/body forces at the
element level) or leaves it to the core.

**C3D8 solid:** does **not** add external loads at the element level;
`F_external` is assembled centrally from concentrated `*Cload` by
`core/loads.external_force`. Distributed/body loads are parsed-not-assembled (a
documented limitation).

## 10. Sign convention

Declared by the class attribute `sign_convention`. `"residual"` means the returned
vector is the `+F_internal` contribution to `R = F_int - F_ext` (the standard
here). `"rhs"` would mean an Abaqus-UEL `RHS = -R`, which must be documented and
converted.

**C3D8 solid, Mode 1, Mode 2:** all set `sign_convention = "residual"` and return
`+F_internal`. The assembler adds `-F_external` itself
(`assembler.py`: `R = R - F_ext`).

**Mode 3 (`uel_adapter.py`):** Abaqus reports `RHS = -R` and `AMATRX = dR/du`.
The adapter converts with the flag `rhs_is_negative_residual` (default `True`):

```
element_residual = -RHS = +R      element_tangent = AMATRX = dR/du
```

so it normalizes to the framework's `+R` and therefore also advertises
`sign_convention = "residual"` (not `"rhs"`) — downstream code sees a standard
residual formulation. The `False` case passes `RHS = +R` through unchanged.

## 11. Verification tests required before the backend is trusted

The levels of the ladder in `docs/verification_strategy.md` that apply, declared
in the class attribute `verification_levels`, plus which are currently passing.

**Per-backend `verification_levels` (from the source):**

| Backend | file | `verification_levels` |
|---|---|---|
| `solid_c3d8_finite_strain` | `formulations/solid_c3d8_finite_strain.py` | `(0,1,2,3,4,5,6,7)` |
| `solid_c3d8_small_strain`  | `formulations/solid_c3d8_small_strain.py`  | `(0,1,2,3,4,5)` |
| `stress_driven_c3d8`       | `formulations/stress_driven_adapter.py`    | `(0,1,3,5)` |
| `uel_direct`               | `formulations/uel_adapter.py`              | `(1,4,5,7)` |

**C3D8 solid, currently passing offline (Abaqus-free):** Level 0 (parser/model
inventory), Level 1 (zero-field residual, `core/verification.zero_field_residual`),
Level 3 (divergence-theorem + linear-stress patch tests,
`tests/cp_c3d8_umat/tangent_fd_check/run_checks.py`), Level 4 (finite-difference
tangent, `core/verification.finite_difference_tangent`; small-strain machine
precision). Levels 5/6/7 (solver comparison / material replay / full residual
replay against Abaqus) are built and documented but **await an Abaqus run** (and,
for the real CP UMAT replay, Intel ifort). See `docs/verification_strategy.md` for
what each level proves, whether it needs Abaqus, and its code realization.

---

**The framework must not claim to support a formulation unless these 11 items are
defined.**
