# Verification Zoo — Manifest Summary

One line per case. Full cards under the category folders; provenance/licenses in
[`sources/`](../../sources). No external code is vendored or executed.

## UMAT cases (`umat_cases/`)

| Case | License | Physics | Backend / Mode | Class | Status |
|---|---|---|---|---|---|
| viscoelasticity | MIT | finite viscoelasticity | `solid_c3d8_finite_strain` · M2 / M1 | B | planned |
| j2_plasticity | MIT | J2 von-Mises plasticity | `solid_c3d8_*` · M2 | B | planned |
| mazars_damage | BSD-3 ⚠ | concrete damage | `solid_c3d8_small_strain` · M2 | B | planned (hold) |
| neohookean_hyperelastic | GPL-3.0 | Neo-Hookean hyperelastic | `solid_c3d8_finite_strain` · M2 | B | reference-only |

## UEL cases (`uel_cases/`)

| Case | License | Physics | Backend / Mode | Class | Status |
|---|---|---|---|---|---|
| uel_elasticity | BSD-3 | linear elasticity (2D/3D Hex8/20) | `uel_direct` · M3 | D | adapter-skeleton |
| uel_hyperelasticity | BSD-3 | finite-strain hyperelastic | `uel_direct` · M3 | D | adapter-skeleton |
| uel_cosserat | MIT | 2D Cosserat/classical | `uel_direct` · M3 (2D) | D | planned |

## Stress-driven cases (`stress_driven_cases/`)

| Case | License | Physics | Backend / Mode | Class | Status |
|---|---|---|---|---|---|
| cp_c3d8_internal | MIT | crystal plasticity (one backend) | `stress_driven_c3d8` · M1 | F | **implemented/offline-verified** |
| external_umat_template | per source | any stress-producing law | `stress_driven_c3d8` · M1 | F | planned (blocked on export) |

## Unsupported cases (`unsupported_cases/`)

| Case | License | Physics | Backend / Mode | Class | Status |
|---|---|---|---|---|---|
| bilinear_czm_cohesive | GPL-3.0 | cohesive (bilinear TSL) | none · M2 | G | unsupported / reference-only |
| uel_hydrogel_coupled | custom | coupled chemo-mechanics | coupled-field · M5 | E | reference-only |
| lsdyna_shell_reference | unverified | user shell (LS-Dyna) | `shell_placeholder` | H | reference-only (wrong solver) |

## Totals

- 4 UMAT + 3 UEL + 2 stress-driven + 3 unsupported = **12 cards**.
- Verified offline **today**: `cp_c3d8_internal` (Mode 1) — plus the framework's
  own truss/beam/mixed/solid suite.
- Everything else is honestly `planned`, `adapter-skeleton`, `unsupported`, or
  `reference-only`, each with a stated **next missing item**.
