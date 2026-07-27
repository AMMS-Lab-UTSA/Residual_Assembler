# Support Matrix

Single source of truth for what the framework supports **today**, spanning both
**internal backends** and the **external verification-zoo candidates**. Honest
statuses only; "offline-verified" means reproduced without Abaqus.

Status vocabulary: `implemented/offline-verified` · `implemented/Abaqus-pending`
· `adapter-skeleton` · `planned` · `unsupported` · `reference-only/license-limited`.

## Internal backends (the framework itself)

| Case | Source | License | Type | Physics | Element / Formulation | Backend required | Current status | Offline verification | Abaqus verification | Notes |
|---|---|---|---|---|---|---|---|---|---|---|
| Stress-driven solid | internal | MIT | Mode 1 | any (field supplied) | C3D8 | `stress_driven_c3d8` | implemented/offline-verified | ✅ `test_assembler.py` (bitwise) | pending ODB (non-uniform state) | Assembly isolated from material |
| Elastic material replay | internal | MIT | Mode 2 | linear elasticity | C3D8 | `solid_c3d8_small_strain` + `isotropic_elastic` | implemented/offline-verified | ✅ kernel match + FD tangent | n/a | Fully runnable offline |
| Crystal plasticity | internal (ngrilli) | MIT | Mode 2 | crystal plasticity | C3D8 | `solid_c3d8_finite_strain` + `crystal_plasticity` | implemented/Abaqus-pending | ✅ kernel/tangent (mock UMAT) | pending ifort+Abaqus | **One** example backend |
| Truss / bar | internal | MIT | Mode 4 | linear axial | T3D2/T2D2 | `truss2` | implemented/offline-verified | ✅ `EA/L` + FD tangent | n/a | Proof of agnosticism |
| Beam / frame | internal | MIT | Mode 4 | Euler-Bernoulli | B31/B33 | `beam2` | implemented/offline-verified | ✅ `PL³/3EI` + FD tangent | n/a | Rotational DOFs |
| UEL direct | internal | MIT | Mode 3 | any (user element) | U1 | `uel_direct` | adapter-skeleton | ✅ sign/FD self-test | pending real UEL | Skeleton adapter |
| Shell | internal | MIT | Mode 4 | shell | S3/S4/... | `shell_placeholder` | unsupported (contract-only) | n/a | n/a | Contract declared, not runnable |

## External candidates (verification zoo)

| Case | Source | License | Type | Physics | Element | Backend required | Current status | Offline verification | Abaqus verification | Notes |
|---|---|---|---|---|---|---|---|---|---|---|
| viscoelasticity | thealanjason | MIT | UMAT | viscoelasticity | C3D8 | `solid_c3d8_finite_strain` (M2/M1) | planned | on export (M1) | required | history-dependent |
| j2_plasticity | UMAT4COMSOL | MIT | UMAT | J2 plasticity | solids | `solid_c3d8_*` (M2) | planned | on export (M1) | required | **non-CP** plasticity |
| mazars_damage | marioruiarruda | BSD-3 ⚠ | UMAT | damage | solids | `solid_c3d8_small_strain` (M2) | planned (hold) | on export (M1) | required | extra permissions file |
| neohookean_hyperelastic | Sina-Taghizadeh | GPL-3.0 | UMAT | hyperelastic | solids | `solid_c3d8_finite_strain` (M2) | reference-only/license-limited | code not reusable | required | re-derive if needed |
| bilinear_czm | harshaa765 | GPL-3.0 | UMAT | cohesive | 2D cohesive | none (not impl.) | unsupported + reference-only | no | required | needs cohesive backend |
| uel_elasticity | bibekananda-datta | BSD-3 | UEL | linear elastic | 2D/3D Hex8/20 | `uel_direct` (M3) | adapter-skeleton | adapter self-test | required | mirrored permissive |
| uel_hyperelasticity | bibekananda-datta | BSD-3 | UEL | hyperelastic | 3D Hex8/20 | `uel_direct` (M3) | adapter-skeleton | adapter self-test | required | unsymmetric AMATRX |
| uel_cosserat | jgomezc1 | MIT | UEL | Cosserat 2D | Quad8/9 | `uel_direct` (M3, 2D) | planned | — | required | micro-rotation DOF |
| uel_hydrogel_coupled | bibekananda-datta | custom | UEL | coupled chemo-mech | coupled (u+μ) | coupled-field (M5) | reference-only | no | required | Mode-5 not implemented |
| uel_generic_demo | CAE Assistant | MIT (demo) | UEL | generic | — | `uel_direct` (concept) | reference-only | no | n/a | incomplete teaser |
| lsdyna_shell_reference | jfriedlein | unverified | USLD/USHL | shell | shell | `shell_placeholder` | reference-only | no | n/a | wrong solver |

## Coverage gaps (documented, not hidden)

| Category | Status | Coverage / plan |
|---|---|---|
| Truss/bar UEL (Abaqus) | no permissive public example found | native `truss2` (Mode 4) covers it |
| Beam/frame UEL (Abaqus) | no permissive public example found | native `beam2` (Mode 4) covers it |
| Cohesive backend | not implemented | GPL CZM reference exists (Mode 2, needs independent backend) |
| Coupled-field (Mode 5) | not implemented | hydrogel UEL reference (custom license) |
| Shell backend | contract only | no permissive Abaqus shell UEL found; LS-Dyna reference only |
| VUMAT (Class C) | none catalogued | explicit-dynamics material path not in scope |

## One-line summary

Offline-verified **today**: stress-driven C3D8, elastic replay, truss, beam,
mixed dispatch (+ the whole `tests/framework` suite). Everything else is honestly
`planned` / `adapter-skeleton` / `unsupported` / `reference-only`, each with a
concrete next step. Crystal plasticity is **one** row, not the table.
