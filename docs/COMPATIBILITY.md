# Compatibility between Residual_Assembler and UMAT-OTI

The two repositories are **separate products**, versioned and released
independently, connected by a versioned contract. Neither vendors the other.

## Current pin

| Component | Version / commit |
|---|---|
| `residual-assembler` | 0.1.0 |
| `umat-oti` | `5fb7b434fc12d0b2cc0ebb101b70cdd9754fbcce` (extra `bridge`) |
| Contract schema | `resasm_umat_transform_v2` |
| Material driver ABI | stdin property vector + strain path; stdout stress and state per increment |

The `bridge` extra pins UMAT-OTI to a **commit**, not a branch. A branch pin
would make this repository's CI result depend on when it ran rather than on what
it contains, and a bridge test that fails only on Tuesdays teaches nobody
anything.

```bash
pip install -e ".[yaml,bridge]"
```

Without the extra, the bridge tests skip and say why. CI installs it, so the
cross-repository contract is genuinely exercised on every push rather than
assumed.

## What crosses the boundary

| Artefact | Producer | Consumer | Validated against |
|---|---|---|---|
| Transformed UMAT object | UMAT-OTI | Residual_Assembler | compilation and primal parity |
| Derivative contract JSON | UMAT-OTI | Residual_Assembler | `resasm_umat_transform_v2` |
| Material-point replay | either | either | stress and state agreement per increment |

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
