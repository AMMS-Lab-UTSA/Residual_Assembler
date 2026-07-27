# UEL Validation Plan

How to validate an external Abaqus **UEL** against the framework's `uel_direct`
(Mode 3) backend. The goal is to confirm that an element's residual/tangent, its
sign convention, and its data layout match what the framework expects — before
any UEL is trusted or wired in.

Scripts: [`scripts/compare_uel_rhs.py`](../../scripts/compare_uel_rhs.py),
[`scripts/run_abaqus_validation.py`](../../scripts/run_abaqus_validation.py).
All skip cleanly when Abaqus is absent.

## 1. What must be verified

| Item | Why it matters |
|---|---|
| **RHS sign convention** | Abaqus reports `RHS = -R`; the framework uses `element_residual = -RHS`. A wrong sign silently negates the residual. |
| **AMATRX convention** | `AMATRX = dR/du`; must match the framework's tangent (and a finite-difference check of the residual). |
| **DOF ordering** | Node-major layout `[u1x,u1y,u1z, u2x,...]`; a permutation corrupts the scatter. |
| **State variable layout** | `SVARS` order/size (`*Depvar`); needed for history-dependent elements. |
| **Element connectivity** | Node order defining the element topology. |
| **Load contribution convention** | Whether body/traction loads are inside the UEL or applied via overlay elements. |

## 2. Instrument the UEL (capture, do not modify the physics)

Add print/write statements inside the UEL to dump, for a chosen element and
increment, a JSON record with exactly these fields:

```
COORDS      -> "coords"        (n_node x ndim)
U           -> "u"             (element DOF vector, node-major)
PROPS       -> "props"
SVARS (in)  -> "svars_before"
RHS         -> "rhs"
AMATRX      -> "amatrx"        (ndof x ndof)
SVARS (out) -> "svars_after"
TIME, DTIME -> "time", "dtime"
JELEM       -> "element"
KINC/NPT    -> integration-point data if available
```

For copyleft/unknown-license UELs, **do not copy the source**; instrument a build
you are licensed to run and export only the *numbers*, or re-implement the element
independently.

## 3. Compare against the framework

```bash
# offline: framework adapter sign/FD self-test
python scripts/compare_uel_rhs.py --framework-selftest

# with an instrumented dump:
python scripts/compare_uel_rhs.py --uel-dump uel_dump.json
```

The script reports: `ndof`, AMATRX shape/symmetry, `element_residual = -RHS`, the
detected sign convention, and an SPD hint for the symmetrized tangent. It reminds
you to confirm DOF ordering, connectivity, and SVARS layout against your model.

## 4. Cross-check paths

1. **Direct (Mode 3)**: feed the same `coords, u, props, svars` into an
   independent element routine and compare `RHS`/`AMATRX` element-by-element.
2. **Assembly (Mode 3 → global)**: scatter the element residuals through the core
   assembler and confirm the global free-DOF residual `~ 0` and reactions vs `RF`
   (`scripts/compare_residuals.py`).
3. **Finite-difference tangent**: perturb `u`, recompute `RHS`, and check
   `-dRHS/du ≈ AMATRX` (the framework's `uel_adapter` already does this for its
   own spring self-test).

## 5. Candidate coverage (see `tests/verification_zoo/uel_cases/`)

| Candidate | License | What it exercises |
|---|---|---|
| Abaqus-UEL-Elasticity | BSD-3 | baseline 2D/3D linear UEL, symmetric AMATRX |
| Abaqus-UEL-Hyperelasticity | BSD-3 | finite strain, **unsymmetric** AMATRX (F-bar) |
| ABAQUS-US (Cosserat) | MIT | extra **micro-rotation DOF** (heterogeneous DOFs) |
| Abaqus-UEL-Hydrogel | custom (ref-only) | **coupled** DOFs (u + μ) — Mode 5, not implemented |

## 6. Acceptance for a UEL to move from `adapter-skeleton` → verified

- sign convention confirmed (`element_residual = -RHS`),
- DOF ordering + connectivity confirmed,
- element `RHS`/`AMATRX` match an independent routine to tolerance,
- global free-DOF residual `~ 0` and reactions match `RF`,
- (if stateful) SVARS layout confirmed and history replayed in order.
