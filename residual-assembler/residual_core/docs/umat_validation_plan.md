# UMAT Validation Plan

How to validate an external Abaqus **UMAT** with the framework. A UMAT answers the
constitutive question only (STRESS, DDSDDE, STATEV at a point); the element
residual is assembled separately. Two complementary checks:

- **Mode 1 (stress-driven)** — assemble the residual from the UMAT's *exported*
  stress field, independent of the material update. This isolates and verifies the
  **assembly**.
- **Mode 2 (material-replay)** — recompute STRESS/STATEV outside Abaqus and
  compare to the ODB, then assemble. This verifies the **material update**.

Scripts: [`scripts/extract_odb_fields.py`](../../scripts/extract_odb_fields.py),
[`scripts/compare_residuals.py`](../../scripts/compare_residuals.py),
[`scripts/compare_umat_replay.py`](../../scripts/compare_umat_replay.py). All skip
cleanly without Abaqus.

## 1. What must be verified

| Item | Why it matters |
|---|---|
| **STRESS** | The Cauchy (or PK) stress in the declared measure/Voigt order. |
| **DDSDDE** | Consistent material tangent; drives Newton convergence and the assembled `K`. |
| **STATEV evolution** | History variables; must evolve increment by increment. |
| **Material parameter mapping** | PROPS order/units → the material's parameters. |
| **Time-increment dependence** | Rate-dependent laws depend on `dtime`. |
| **Finite-strain conventions** | Stress measure (Cauchy/PK2), kinematic input (F), objective rates. |

## 2. Mode 1 — stress-driven (verifies assembly; offline-capable)

```bash
# 1. run the job (needs Abaqus + ifort for the UMAT)
python scripts/run_abaqus_validation.py --job case --inp model.inp --user umat.for
# 2. export the integration-point stress (inside Abaqus python)
abaqus python scripts/extract_odb_fields.py --odb case.odb --out fields.json
# 3. assemble the residual OUTSIDE Abaqus and compare (no Abaqus needed here)
python scripts/compare_residuals.py --model model.json --fields fields.json
```

Expected: free-DOF residual `~ 0` at equilibrium; reactions match `RF`
(sign-aware). This is exactly the path the internal C3D8 crystal-plasticity case
already passes offline (`tests/verification_zoo/stress_driven_cases/cp_c3d8_internal`).
**The material model is irrelevant to this check** — that is the agnostic claim.

## 3. Mode 2 — material replay (verifies the update)

### 3.1 Single material point / single element
Drive the material through the **same** kinematic history the ODB saw, and compare
STRESS/STATEV per increment:

```bash
abaqus python scripts/extract_odb_fields.py --odb case.odb --out ref.json   # per-frame history
python scripts/compare_umat_replay.py --reference ref.json                  # compiled UMAT (needs Abaqus/ifort)
python scripts/compare_umat_replay.py --reference ref.json --python-material mymod:MyMaterial  # offline, independent code
```

### 3.2 History-dependent materials — the hard rule
> For plasticity, viscoelasticity, damage, and crystal plasticity, a
> **final-step-only replay is invalid**. The load history must be replayed
> increment by increment, propagating STATEV. `compare_umat_replay.py` iterates
> the exported frame sequence and the framework's `StateManager` enforces the
> committed/trial discipline.

Path-independent laws (linear elastic, hyperelastic) are exempt: a single step
suffices.

## 4. Licensing constraint on replay

- Permissive UMAT (MIT/BSD): may be compiled and replayed (`--python-material`
  only if you *independently* re-implemented it).
- Copyleft / unknown UMAT: **do not copy or link**. Use Mode 1 (export stress,
  assemble) which needs only the *numbers*, or re-derive the law from published
  equations for an independent Mode-2 replay.

## 5. Candidate coverage (see `tests/verification_zoo/umat_cases/`)

| Candidate | License | Path | History? |
|---|---|---|---|
| viscoelasticity (MIT) | permissive | M1 + M2 | yes |
| J2 plasticity (MIT) | permissive | M1 + M2 | yes |
| Mazars damage (BSD-3 ⚠) | permissive* | M1 + M2 | yes (softening) |
| Neo-Hookean (GPL) | reference-only | M1 only (code not reusable) | no |
| crystal plasticity (MIT, internal) | permissive | **M1 offline-verified** | yes |

## 6. Acceptance for a UMAT case to move from `planned` → verified

- **Mode 1**: assembled free-DOF residual `~ 0`, reactions match `RF` — *offline*
  once the stress export exists.
- **Mode 2**: replayed STRESS/STATEV match the ODB per increment to tolerance,
  with the full history replayed in order.
