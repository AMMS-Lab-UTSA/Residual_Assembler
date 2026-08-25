# Changelog

Notable changes to Residual_Assembler. This project follows
[semantic versioning](https://semver.org/); while the major version is 0 the
public interfaces may still change.

## [Unreleased]

### Added
- Public project files: `AUTHORS.md`, `SECURITY.md`, `CODE_OF_CONDUCT.md`,
  `CONTRIBUTING.md`, `CITATION.cff`, `codemeta.json`, `.zenodo.json`,
  `THIRD_PARTY_LICENSES.md`, `.gitattributes`, this changelog.
- `tools/audit_repository_standards.py`, run in CI: required files, no tracked
  build products, no secrets, no absolute home paths.
- Explicit pytest markers (`unit`, `integration`, `slow`, `network`, `fortran`,
  `abaqus`, `arc`, `publication`) so a missing Abaqus installation cannot
  silently shrink the offline suite.

### Known gaps
- **No licence is declared for this repository's own code.** Default copyright
  therefore applies. See `AUTHORS.md`; the audit reports this until it is fixed.
- The sensitivity beat does not yet run through the C3D8 assembly path. See
  `STATUS.md`.

## [0.1.0]

### Added
- Fresh-clone reproducibility: external Abaqus sources are pinned git submodules
  organised by licence tier, never vendored. `sources/SUBMODULES.md` records the
  policy and every pinned commit.
- `scripts/init_permissive_sources.sh` fetches only the permissive tier and
  refuses any path outside it.
- `scripts/verify_source_submodules.py` and `scripts/verify_clean_clone.sh`.
- `REQUIRE_EXTERNAL_TEST_SOURCES=1` turns skips caused by absent external
  sources into failures, so a misconfigured checkout goes red rather than
  green-but-smaller.
- C3D8 residual assembly verified from exported Abaqus ingredients.
