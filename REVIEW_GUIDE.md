# Review Guide

A guided tour of the **model-agnostic residual assembly framework** for a
reviewer. Read top to bottom; every offline step runs on Python 3 + numpy with no
Abaqus.

## 1. What the framework is

It assembles the global finite-element **residual**

```
R(y, q, a, t) = F_internal(y, q, a, t) - F_external(t) + F_constraints(y, t)
```

*outside* the FE solver, for arbitrary element formulations and materials. The
core assembler contains **no** element or material mathematics: it dispatches each
element to a registered **formulation backend** and each material point to a
registered **material backend**, then scatters. See
[`residual_core/docs/architecture.md`](residual_core/docs/architecture.md).

## 2. Why it is model-agnostic (not just asserted)

Argued from the code in
[`residual_core/docs/minimality_proof.md`](residual_core/docs/minimality_proof.md):
no physics in the core, heterogeneous per-node DOFs, backends that *declare* their
capabilities, and a requirements engine that asks for the **minimum** missing
input. The same assembler drives stress-driven solids, elastic replay, crystal
plasticity, truss, beam, and a mixed model unchanged.

## 3. Implemented backends

Seven formulation backends and three materials (`resasm backends`):

- `stress_driven_c3d8` (Mode 1), `solid_c3d8_small_strain` / `solid_c3d8_finite_strain`
  (Mode 2), `truss2` / `beam2` (Mode 4), `uel_direct` (Mode 3, skeleton),
  `shell_placeholder` (contract only);
- materials: `isotropic_elastic`, `umat`, `crystal_plasticity`.

**Crystal plasticity is one backend, not the purpose of the framework.** Full
status: [`residual_core/docs/support_matrix.md`](residual_core/docs/support_matrix.md).

## 4. External UMAT/UEL cases (verification zoo)

Public Abaqus user subroutines were catalogued across diverse physics —
**viscoelasticity, J2 plasticity, damage, hyperelasticity, cohesive, coupled
chemo-mechanics** — plus crystal plasticity. **No external code is vendored**;
copyleft/unknown/different-solver sources are reference-only.

- Provenance + licenses: [`sources/external_manifest.md`](sources/external_manifest.md),
  [`sources/external_manifest.json`](sources/external_manifest.json),
  [`sources/LICENSE_AUDIT.md`](sources/LICENSE_AUDIT.md).
- Verification cards: [`tests/verification_zoo/`](tests/verification_zoo/) — one
  `CASE.md` per case with source, license, backend, required inputs, expected
  outputs, offline vs. Abaqus tests, status, and the **next missing item**.

## 5. How they are classified

Each case has a **residual class** (A–H) and a **minimum verification mode** (1
stress-driven / 2 material-replay / 3 direct-UEL / 4 native formulation) — see
[`tests/verification_zoo/README.md`](tests/verification_zoo/README.md) and
[`manifest_summary.md`](tests/verification_zoo/manifest_summary.md).

## 6. What is verified offline (today)

```bash
# one-shot guided demo
bash scripts/demo_framework.sh

# or the individual suites
python tests/framework/test_assembler.py               # core == kernel (bitwise)
python tests/framework/test_truss2_backend.py          # EA/L
python tests/framework/test_beam2_backend.py           # PL^3/3EI
python tests/framework/test_mixed_model_dispatch.py    # physics-blind dispatch
python tests/framework/test_dof_manager_mixed.py       # heterogeneous DOFs
python tests/framework/test_requirements_negative.py   # minimum missing input
python tests/framework/test_public_api.py              # facade only
python tests/framework/test_neutral_io.py              # neutral JSON round-trip
python -m residual_core.ui.examples                    # runnable API examples
```

Transcript of real output:
[`residual_core/docs/demo_transcript.md`](residual_core/docs/demo_transcript.md).

## 7. What still needs Abaqus

Everything that needs a real ODB or a compiled user subroutine is **ready-to-run
and clearly marked pending**, never claimed as passing. The order to close each
gap is in
[`residual_core/docs/abaqus_validation_roadmap.md`](residual_core/docs/abaqus_validation_roadmap.md).
Ready-to-run scripts (skip cleanly without Abaqus):

```bash
python scripts/run_abaqus_validation.py --job J --inp model.inp --user umat.for
abaqus python scripts/extract_odb_fields.py --odb J.odb --out fields.json
python scripts/compare_residuals.py   --model model.json --fields fields.json
python scripts/compare_umat_replay.py --reference ref.json
python scripts/compare_uel_rhs.py     --uel-dump uel.json
```

Plans: [`umat_validation_plan.md`](residual_core/docs/umat_validation_plan.md),
[`uel_validation_plan.md`](residual_core/docs/uel_validation_plan.md).

## 8. How to run demos/tests

- **Demo**: `bash scripts/demo_framework.sh` (offline).
- **API**: `from residual_core import ResidualProblem` — see
  [`residual_core/docs/user_interface.md`](residual_core/docs/user_interface.md).
- **CLI**: `python -m residual_core.ui.cli --help` (or `resasm` if installed).
- **Tests**: the `tests/framework/*.py` scripts above; each prints `PASS`/`FAIL`.

## 9. Guardrails honored

OTILib is the real HYPAD backend (external dependency; reported cleanly when
absent, no silent fallback). Dual1 is only a legacy first-order smoke test. No
fake OTI is bundled and no 150-parameter model, unsupported physics, or
restricted/unknown-license code is integrated. Crystal plasticity remains one
verified example backend among many.
