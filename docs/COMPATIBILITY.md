# Compatibility between Residual_Assembler and UMAT-OTI

The two repositories are **separate products**, versioned and released
independently, connected by a versioned contract. Neither vendors the other.

## Current pin

| Component | Version / commit |
|---|---|
| `residual-assembler` | 0.1.0 |
| `umat-oti` | `4822108cb977160f5b66d46d50753bb9ad41e848` (extra `bridge`) |
| Shared contract | version 3.0.0; `schemas/contract_lock.json` equals UMAT-OTI's byte for byte |
| Transform generation | `da1f183708c19072` (`schemas/transform_generation.json`) |
| Compiled provider ABI | `UMAT`, `UMAT_OTI_EVAL`, `UMAT_OTI_MARCH`, `UMAT_OTI_EVAL_TOTAL`, described by the completed contract (`Mapping.json`) |
| Contract schema (driver path) | `resasm_umat_transform_v2` |
| Material driver ABI | stdin property vector + strain path; stdout stress and state per increment |

The `bridge` extra pins UMAT-OTI to a **commit**, not a branch. A branch pin
would make this repository's CI result depend on when it ran rather than on what
it contains, and a bridge test that fails only on Tuesdays teaches nobody
anything.

```bash
pip install -e ".[yaml,bridge]"
```

CI checks out UMAT-OTI at the same commit beside this repository and installs
it from there, because the connected tests also read its source tree
(contracts and UMAT sources).

Without the extra, the bridge tests skip and say why. CI installs it, so the
cross-repository contract is genuinely exercised on every push rather than
assumed.

## What crosses the boundary

| Artefact | Producer | Consumer | Validated against |
|---|---|---|---|
| Transformed UMAT object | UMAT-OTI | Residual_Assembler | compilation and primal parity |
| Derivative contract JSON | UMAT-OTI | Residual_Assembler | `resasm_umat_transform_v2` |
| Material-point replay | either | either | stress and state agreement per increment |

## What this repository needs and does not yet get

Three things cross the boundary today as numbers and one does not cross at all.
Named here so that neither side has to guess.

### From UMAT-OTI

1. **The transform fingerprint, on every artefact.** Already carried by
   `export_residual_fixture.py` as `transform_fingerprint`, and now READ:
   `verified_fixture.load()` refuses a fixture that does not match the current
   store. The request in the other direction is that the promoted collection
   under `umat/` carry it in the same field and at the same value — it
   currently records `e4257779bd847cc4` in `contract.json` while the fixtures
   from `pass11` carry `b0d27ee53c630500`, so the two artefacts of one run
   disagree about which build produced them.
2. **`dsigma/dq` and `dstatev/dq`, per integration point.** The assembler can
   now assemble `dR/dq` (`core/state_sensitivity.py`) and has nothing verified
   to put in it: the `state_sensitivity` layer of the attribution reports
   `not_established` for every committed fixture. The shape needed is
   `(n_ip, ntens, n_state)` per element, in the same Voigt order and the same
   integration-point order as `STRESS`, at a NAMED increment — a state
   derivative without the state it was taken at is not usable.
3. **Which reading of DDSDDE a converted build's tangent is in.** Measured
   here per fixture (`tangent_convention`) because the fixture does not say,
   and the answer is a property of the routine rather than of the deck's
   `NLGEOM`. If the transformation already knows — it has the source — saying
   so in the fixture would turn a measurement into a fact and would let a
   material whose window is isochoric be assembled at all.

### From the GUI

1. **`Diagnosis.as_dict()` is the interface, not the prose.** Keys: `ok`,
   `complete`, `blame`, `not_established`, `findings`; each finding carries
   `layer`, `status` (`holds` / `fails` / `not_established`), `detail`,
   `measured` and `would_establish`. `ok` means nothing that was CHECKED
   failed; `complete` means every layer was checked and held. A GUI that
   renders `ok` as a green tick is reporting unchecked layers as passing.
2. **Show the sweep, not the best number.** `Sweep.as_dict()` carries
   `steps`, `errors`, `plateau`, `plateau_span`, `flat` and
   `worst_components`. A single error is not evidence, and an error that is
   flat in the step size is a wrong formula rather than a truncation error —
   `flat: true` is the one field that must never be collapsed into a
   tolerance badge.
3. **A layer with no evidence must not be drawn as a layer that passed.**
   `not_established` needs its own colour.

## Bumping the pin

1. Confirm the UMAT-OTI commit is pushed and its CI is green.
2. Update the `bridge` extra in `pyproject.toml` to the new commit.
3. Run `python -m pytest -q -m "not abaqus and not arc and not network"`.
4. Update the table above and note the change in `CHANGELOG.md`.

Bump deliberately. An automatic bump would remove the only place where an
incompatible contract change is noticed before a user meets it.

## Deprecation policy

While either package is at major version 0, contract fields may be added freely
and removed only after they have been reported as deprecated in at least one
release. A schema change that removes or retypes a field is a breaking change
and requires a major-version bump of the schema name itself
(`resasm_umat_transform_v3`), not a silent redefinition of v2.

## A pin can be invalidated by a history rewrite

Pinning a commit assumes that commit stays reachable. Rewriting the pinned
repository's history -- to correct commit authorship, for instance -- replaces
every SHA on the branch, and the old one stops existing on the remote. `pip`
then fails to resolve the dependency, which is loud and immediate rather than
subtly wrong, but it does mean the pin has to be updated in the same batch as
any such rewrite.

If you rewrite history in UMAT-OTI, bump the `bridge` pin here before pushing.
