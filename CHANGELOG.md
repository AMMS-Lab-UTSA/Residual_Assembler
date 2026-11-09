# Changelog

Notable changes to Residual_Assembler. This project follows
[semantic versioning](https://semver.org/); while the major version is 0 the
public interfaces may still change.

## [Unreleased]

## [0.2.0] - 2026-11-09

### UMAT-OTI 1.2.1 (2026-11-06)
- **Companion pin.** The `bridge` extra (and therefore CI) uses UMAT-OTI 1.2.1,
  tag `v1.2.1`, commit `bb099d85b080a295a443773c823305671164653f`; it was
  `95b27ff`.
  The history of both repositories was rewritten on 2026-11-07; see
  [docs/HISTORY_REWRITE_2026-11.md](docs/HISTORY_REWRITE_2026-11.md).
- **Contract files.** `schemas/transform_generation.json` and
  `schemas/contract_lock.json` are copied byte for byte from UMAT-OTI 1.2.1:
  transform generation `80a2015b780b1f25` (harness `56e3725136cb897e`).
  Contract 5.0.0 is unchanged.
- **The two verified fixtures re-frozen in Abaqus at `80a2015b780b1f25`.**
  `scripts/regenerate_recovery_fixture.py` re-ran both against UMAT-OTI 1.2.1:
  the transform, the original and transformed Abaqus 2021 jobs and the
  Jacobian-matched control, three jobs per fixture, one after another (six in
  all). The elasticity experiment comes from its retained pass25 row; J2 has
  no row of its own and is templated from the elasticity row's manifest with
  the archived J2 material. Both decks are byte-identical to the archived
  ones (docs/evidence/recovery_fixtures.md, 2026-11-06 section). Both are
  verified (4 of 4 chosen states judged entrywise, every entry agrees). Every
  number is identical to the fixtures frozen at `c61511ea0ffe58e2`; only
  `generated` and `transform_fingerprint` changed. The residual check is
  unchanged: elasticity 1 held, 2 not established; J2 32 held, 3 not
  established (`docs/evidence/final_refreeze_fixture_*_check.json`).
- **streamlit >= 1.61** for the `gui` extra: the GUI passes `width="stretch"`
  to data frames, which older versions reject (the floor said 1.30), as in
  UMAT-OTI 1.2.1.
- **Benchmarks without the book code.** The seven claim-4 benchmark contracts
  whose sources are the Dunne & Petrinic book code, which UMAT-OTI no longer
  ships, end as "source not available"; the Abaqus elastic benchmark skips
  without elastic.f (recorded in docs/evidence/integrity_review.json).
- **Docs.** COMPATIBILITY.md, USAGE_REPORT.md, HANDBOOK.md and
  which_layer_did_it_come_from.md give the new pin and generation.
- **Scrub.** Machine paths in shipped records, examples and expected outputs
  are repository-relative or `<home>`; agent names in comments, docstrings
  and record paths name components or roles; the agent-process files
  `.councilignore` and `docs/LOCAL_COUNCIL.md` are removed.

### The deck Abaqus ran
- `resasm history` (and `resasm request` through it) reads a deck with
  `*User Material, ..., UNSYMM` and with results-file or printed-output
  requests (`*File Format`, `*Node File`, `*El File`, `*Energy File`,
  `*Node Print`, `*El Print`, `*Energy Print`) instead of refusing it. Neither
  enters the residual: UNSYMM is how Abaqus stores its own Newton matrix (the
  replay assembles the provider's DDSDDE as returned), and the output requests
  choose what Abaqus writes. Both are listed in `run_report.txt` under "Deck
  lines accepted and not applied". Such a deck reads into exactly the model
  it reads into without them (`tests/replay_history/test_deck_output_requests_and_unsymm.py`).

### Connected workflow (2026-09-18)
- **`resasm history`**: total-history sensitivities of small-strain C3D8
  analyses for any UMAT-OTI provider (entry point `UMAT_OTI_EVAL_TOTAL`),
  with prescribed displacements, many increments, sparse assembly, von Mises
  and volume-weighted outputs, weighted shares, `--reequilibrate` and
  `--verify tangent|fd`. `resasm request` hands every readable model outside
  its bounded scope to it and says so. Both presentation cantilevers (J2,
  1,536 C3D8 / 40 increments; FCC, 384 C3D8 / 25 increments) run through it;
  see `docs/REPLAY_HISTORY.md` and `examples/cantilevers/`.
- **Independent checks** of those sensitivities: whole-model central
  differences of the ORIGINAL UMAT, Abaqus reruns, and Euler's homogeneity
  identity at every increment (`tests/replay_history/test_history_example.py`).
- **GUI**: the Sensitivity Request / Solve screen ticks parameters from
  `Mapping.json`, chooses outputs and region, and shows full-field results.
- **Clean-install gate** (`scripts/clean_install_gate.py`) runs on any clean
  branch (`--branch`) and adds the full-size J2 cantilever (`--cantilever`).
- **Presentation claims** reproduced claim by claim (`verification/run_all.py`,
  `docs/VERIFICATION_RECORD.md`), including what did not reproduce.
- Shared contract 5.0.0 (MAJOR, from UMAT-OTI 7cecbb5): terminal state
  `primal_control_not_decided` (INTERNAL) and `primal_decided_by` value
  `jacobian_matched_not_decided`. The routine-level replay agreed and the Abaqus
  control produced no comparison, so `primal_agreed` is null: NOT ESTABLISHED,
  neither a measured disagreement nor a pass. A 4.x producer is now refused at
  the version handshake. The lock (combined digest `82db77bf...`) and the four
  schemas are copied byte for byte; both fixtures keep verifying unchanged.
- Shared contract 4.1.0: `finite_history.evidence.primal_decided_by` is declared
  (optional); 4.0 records are still read, with that field NOT ESTABLISHED.
- Shared transform generation re-frozen at `c61511ea0ffe58e2` (previously
  `50ae446f12de66d8`) for UMAT-OTI corpus pass25, with the two fixtures
  regenerated in Abaqus 2021. Every number is identical; only `generated` and
  the fingerprint changed. The residual check is unchanged (elasticity 1 held /
  2 not established, J2 32 held / 3 not established).
- Shared transform generation re-frozen at `50ae446f12de66d8` (previously
  `830e5ee95ce99cd2`) for UMAT-OTI corpus pass24, with the two fixtures
  regenerated in Abaqus 2021. Every number is identical; only `generated` and
  the fingerprint changed. The residual check is unchanged (elasticity 1 held /
  2 not established, J2 32 held / 3 not established).
- Shared transform generation re-frozen at `830e5ee95ce99cd2` (previously
  `a4f0ea8c9d124f18`, `dbe9f928191e1d43`) with the two current fixtures
  regenerated in Abaqus. The numbers are unchanged; the producer's entrywise
  tangent gate now records the states it judged (`states_agreeing` null), and
  `verified_fixture` describes such a fixture accordingly
  (`docs/evidence/final_refreeze.md`).
- The `bridge` extra pins UMAT-OTI `95b27ff`; CI checks the companion out at
  that commit beside this repository (`docs/COMPATIBILITY.md`).

### Added
- **Fixture provenance gate.** `verified_fixture.load()` refuses a fixture
  whose `transform_fingerprint` is not the current store's
  (`b0d27ee53c630500`). The four fixtures this repository verified against
  were frozen under `ff94800b1884bcc0`, so every comparison made against them
  was a comparison with a transformation two builds old, and nothing said so.
  The set is re-frozen from `corpus_run/pass11`: nine cases, nine
  repositories, both kinematics, a moving state, two rate-dependent
  materials, a plane-stress case the C3D8 kernel refuses, and one case whose
  two builds disagreed.
- **`tangent_convention()`** measures which reading of DDSDDE a fixture's own
  numbers support — the plain `D:Deps` or Abaqus's finite-strain
  `D:Deps - sigma tr(Deps)`. It is a property of the ROUTINE and not of the
  step's `NLGEOM` flag: `irfancn umat_elastic` under `NLGEOM=YES` satisfies
  the plain reading to 3.4e-16 and misses the other by 5.5e-03; the
  neo-Hookean models under the same flag are the other way round.
- **Nine-layer attribution** (`residual_core.diagnose`), replacing the five
  stages: the original UMAT, the OTI transformation, the constitutive
  tangent, the state sensitivity, the tensor mapping, the element
  integration, the residual assembly, the boundary conditions, and the global
  DOF mapping. Every verdict carries its measurement; `not_established` is a
  third answer distinct from holding and says what would settle it. See
  `docs/which_layer_did_it_come_from.md`.
- **`core/state_sensitivity.py`**: `dR/dq` assembly, block-sparse by
  construction — a parameter is shared by every integration point and a state
  is not.
- **`core/deck_replay.py`**: reconstructs the motion a fixture's own deck
  prescribes, which is what makes the boundary and DOF layers checkable. It
  also established that Abaqus's `NLGEOM` `STRAN` is the mid-point
  (Hughes-Winget) accumulation and not `ln V`: measured at 0.9% extension,
  `ln V` misses by 8.3e-08 and the accumulation reproduces it to 2.2e-16.
- **`core/finite_difference.py`**: one plateau implementation for the whole
  repository, including flatness detection — an error that does not move with
  the step size is a wrong formula, and it reads exactly like a converged
  plateau to anything that only asks whether the best number is small.
- Integration tests driven by the verified fixtures: the strain reconstructed
  from a fixture's own deck against the strain Abaqus recorded (2.9e-16 small
  strain, 2.1e-14 finite); global `dR/du` across eight elements over a
  plateau; `dR/dq` against a closed-form state derivative (plateau of 3,
  best 4.3e-10); global `dR/dp` and `du/dp` against re-solving the model
  (plateaus of 5, best 3.3e-14 and 6.1e-11).
- Public project files: `AUTHORS.md`, `SECURITY.md`, `CODE_OF_CONDUCT.md`,
  `CONTRIBUTING.md`, `CITATION.cff`, `codemeta.json`, `.zenodo.json`,
  `THIRD_PARTY_LICENSES.md`, `.gitattributes`, this changelog.
- `tools/audit_repository_standards.py`, run in CI: required files, no tracked
  build products, no secrets, no absolute home paths.
- Explicit pytest markers (`unit`, `integration`, `slow`, `network`, `fortran`,
  `abaqus`, `arc`, `publication`) so a missing Abaqus installation cannot
  silently shrink the offline suite.

### Changed
- **`*STEP` is no longer flattened away.** The parser records which step each
  `*BOUNDARY` and `*CLOAD` was written in and its `OP=NEW`/`OP=MOD`, and
  `constraints.active_boundaries(model, step)` applies Abaqus's rule. The
  verification decks carry four steps each opening with `OP=NEW`; read flat, a
  deck's first step was driven by its fourth step's displacements, silently.
  `dirichlet_dofs`/`partition` take an optional `step` and keep their previous
  behaviour when it is not given.
- Licensed **GPL-3.0-only**, matching the companion product UMAT-OTI and the
  OTILib algebra both build on. The `sources/` licence tiers keep their
  `update = none` policy but their rationale changes: they are a distribution
  and obligation boundary, not a firewall around a non-copyleft core.

### Known gaps
- No J2 elastoplastic fixture exists at the current transform fingerprint: the
  bundled control was last verified at `ff94800b1884bcc0` and `pass11` carries
  no plasticity case. Filed as
  `abaqus_queue/requests/A3_j2_control_refresh.json`. Every committed fixture
  is elastic, hyperelastic, viscous or growth, so the one family whose tangent
  changes between predictor and corrector is unrepresented.
- The finite-strain ELEMENT tangent has still never been compared with
  Abaqus's. The material-point convention is now measured; the element
  consequence needs reaction derivatives from a real run
  (`A3_finite_strain_element_tangent.json`).
- The integration-point ordering is internally consistent and unverified
  against Abaqus's export index; every offline check uses a uniform field,
  which is provably blind to it (`A3_integration_point_ordering.json`).
- `dR/dq` is verified against a closed-form state derivative of the pure-Python
  J2 reference, not against a transformed UMAT's own DSTATEV sensitivities.
  The `state_sensitivity` layer reports `not_established` for every committed
  fixture, because none carries a `dR/dq`.

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
