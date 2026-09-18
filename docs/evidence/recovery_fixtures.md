# Recovery Fixtures: Current Runs And Historical Evidence

Date: 2026-09-18. Full RA offline suite: **419 passed, 0 failed, 11 existing
skips, 0 errors**, 57.05 seconds. No new skips, xfails, weakened tolerances,
mocked mathematics, commits, branches, pushes, subagents, or corpus batch runs.
Existing dirty work was preserved. No producer code, generation, lock, production
fixture guard, or pinned J2 source-hash check was changed in this task.

## Scope And Classification

The earlier result was 325 passed, 56 failed, 11 skipped, 5 errors. The exact
failure list and JUnit in recovery_finite.md were inspected: all 61 failing
entries hit stale generation rejection, including a collection error hiding
the parametrized deck tests. History d3bffa2 was inspected: its legitimate
pass12 regeneration was the precedent, not permission to relabel old numbers.

The current producer remains **6aa20d22e37f14c9**, verified by the actual
transform_fingerprint() function before and after both runs. The ten original
fixtures remain **94a92c01814f107a**, with every byte/hash preserved.

Historical tests explicitly read that recorded generation and assert it is
not current. The test-only historical loader first proves that the production
default loader rejects each archived fixture. No production acceptance uses
this helper. The archived ten-case count, multi-author coverage, plane-stress,
growth, finite-strain, J2 cycle, and failed-primal viscoelastic attribution
example remain covered. No hard example was removed or relabelled.

| Test slice | Evidence it tests |
| --- | --- |
| test_a_fixture_is_held_to_the_rule_that_froze_it | Historical loader claims and deliberately damaged historical payloads |
| test_a_step_is_not_decoration | Recorded four-step J2 deck semantics |
| test_a_verified_deck_drives_the_global_assembly | Historical deck reconstruction, playback, and assembly arithmetic, not a live constitutive solve |
| test_a_verified_umat_drives_the_assembler | Historical tensor mapping, integration, and failure attribution |
| test_the_assembled_residual_derivative_is_differenced | Historical finite-difference assembly checks and original ten-case coverage denominator |
| test_the_assembled_state_derivative_is_differenced | State-derivative attribution using recorded inputs |
| test_which_layer_the_failure_came_from | Historical tangent conventions and known failed-primal subject |
| Contract schema/identity/gate tests | All ten archived contract examples, including legacy schema cases |
| New current fixture/CLI regressions | Two genuinely regenerated current fixtures, default production loader, independent references, unchanged source/deck/count identities |

The current set is deliberately two materials, not a ten-material capability
claim. Elasticity tests the public assembly-consistency workflow. Pinned J2
adds nontrivial state and four material-parameter derivatives against the
separate producer J2 reference. The other eight historical materials have not
been established under current code. Neither collection is the corpus.

## Inputs And Actual Execution

Targeted lookup found all nine external sources in discovery_cache, matching
their archived source SHA-256 values. Bundled J2 was absent at its former cache
path but found at imq-umat-recovery/UMATs/UMATs/generic_ps/j2_props.f with the
exact pinned digest. Thus there is no missing-source blocker for this bounded
repair. No external dependency was presumed unavailable.

Elasticity uses pass12's retained manifest from
corpus_run/_superseded/pass12/results/store_verification.jsonl. Its author's
UMATmodel.inp supplies E=210000, nu=0.3. J2's archived material supplies
E=200000, nu=0.3, SIGY0=250, H=2000, NSTATV=1. The J2 manifest uses the retained
four loading segments, with its material/name restored from the archived
fixture. **generate_deck(manifest) must match the entire archived deck byte for
byte before transformation or solving.** Both checks passed. The four segments
are uniaxial extension, shear, reversal, and a ten-time-unit hold, with
10/10/10/5 increments. These are the original harness experiments; material
constants are author-derived, not an assertion that the loading is an author's
full-sized model.

scripts/regenerate_recovery_fixture.py calls the existing transform_one /
run_transformation API with compilation enabled, writes a new TransformStore,
then calls only verify_one with frozen loading, followed by the documented
export_residual_fixture.freeze. It never calls verify_store_in_abaqus.main or
run_batch. The actual original and current transformed UMATs execute in licensed
Abaqus with ifort. The existing independent original-source FD replay verifies
the OTI tangent. No old numerical verification result is used as a new result.

Exactly four Abaqus analysis jobs ran sequentially, two per selected material:
imqrf_original and imqrf_transformed, in separate material/build directories
under workspace/imq_abaqus/recovery_fixtures. A wrapper changes only job names
and creates same-byte filename aliases for the verifier/freezer's conventional
original/transformed filenames. It does not replace solver, comparator,
mathematics, or evidence gates. Jobs and reports are in each run's jobs.json.
All commands were synchronous; no detached jobs or whole-corpus verification.

Abaqus wrote successful completed analyses for all four runs. The unchanged
job-status classifier also records the known post_analysis_wrapup_failure and
process_exit_code_1 warnings. Success is based on the actual analysis status,
complete finite probe history and independent comparisons, not process code.

## Measured Evidence

| Measurement | Isotropic Elasticity | Bundled J2 |
| --- | ---: | ---: |
| Completed increments, each build | 35 | 35 |
| Integration points per increment | 8 | 8 |
| Full records, each build | 280 | 280 |
| Finite values scanned, each build | 12040 | 12040 |
| Frozen records at named IP 1 | 6 | 35 |
| Frozen steps | Original six-record window | All four steps |
| Worst primal stress/state relative difference | 0 / 0 | 0 / 0 |
| Resolved primal components | 968 | 1312 |
| Unresolved primal components, not counted as independent agreement | 700 | 612 |
| Independent tangent states agreeing / checked | 4 / 4 | 2 / 2 |
| Worst accepted tangent relative error | 1.069e-14 | 8.250e-11 |
| Worst-state accepted FD plateau | 0.01 through 0.1 (2 sizes) | 1e-6 through 1e-4 (3 sizes) |
| Acceptance gates explicitly true | 6 / 6 | 6 / 6 |
| Current RA CLI held / failed / not established | 1 / 0 / 2 | 32 / 0 / 3 |

Both full-history groupings and five-field identities (element, point, step,
increment, time) are retained. Regressions check unchanged full counts, sources,
material data, deck bytes/digest, and frozen lengths against the archive.
Combined current scope: two materials, 560 original plus 560 transformed full
records, 48160 scanned finite values across both sides, six independently
verified tangent states. No unresolved quantity is counted as an agreement.

Elasticity's RA dR/du check is explicitly assembly consistency using its
recorded linearization, not independent constitutive differentiation. Its
dR/dq and dR/dp remain not established. J2 recomputes stress from the previous
state through umat_oti.validation.j2_reference and verifies dR/du, dR/dq,
dR/dE, dR/dNU, dR/dSIGY0, dR/dH. Three comparisons at an elastic increment
have exactly zero state/yield/hardening derivatives and are deliberately not
claimed as passes. The producer's two J2 tangent states are activated loading
states; this does not claim independent producer tangent coverage of every
constitutive branch.

Full current reports: recovery_fixture_elasticity_verification.json and
recovery_fixture_j2_verification.json. RA public reports:
recovery_fixture_elasticity_check.json and recovery_fixture_j2_check.json.
The corresponding *_hashes.json files retain every generated artifact digest,
including transformed source, compiled support, probes, replay results, and
ODB. All 129 elasticity and 111 J2 artifact hashes were rechecked successfully
after the full suite. Source code is not redistributed into these evidence reports.

| Artifact | SHA-256 |
| --- | --- |
| Current elasticity fixture | f9ebddcf4bb319319f447b6fd1d9448e9cbb4d98650fe3ab12362235040a6124 |
| Current J2 fixture | 881479cf8e1555c6635572b616e86734433bb97cc3e7c727432c09b51b032447 |
| Elasticity original source | 27cda337a45be53431149436fbad42bd19da0bb5ec2a68f6370c17b608326a86 |
| J2 original source | 4362e759c4642cb2f74e20c50b86a9472b1e80d6a7fef663bf8faa3dfa8de1c8 |
| Elasticity deck | 148a9ba1eb9cc250940ec7cca532602a66a8bad1bfb363e8933a7816f0018d92 |
| J2 deck | 5bca9c51a7f5aab24d31950ac17a0e70be7d96f65669193cdc76e95d50fe332a |
| Elasticity verification report | d1ad3c3944a7d32660b86bc7ef5a245d5371815d1957b09fc89b8c91aec6cd89 |
| J2 verification report | c9636ffed9175bfce3169bccbb4de502f02575d53b96bbdad8c6c4d4c1f10172 |

## Historical Archive Map

The original recovery_evidence_inventory.json is an immutable historical
snapshot: its 767-file scope and original paths/hashes have NOT been rewritten
to describe this later state. This table is its relocation addendum for RA.
Every old tests/fixtures/verified/NAME now exists byte-identically at
tests/fixtures/historical/94a92c01814f107a/NAME. The archive-hash regression checks
all ten against the original inventory. Current replacements exist only for
elasticity and J2 at their former operational paths; the other eight are
historical-only, not deleted evidence. UMAT inventory entries remain untouched.

| NAME | Original SHA-256, unchanged in archive |
| --- | --- |
| compresibleneohookean--d78e61195f.json | ba9e637e339c37aca90c1aa63fc8e9e72cd21f9326a07336249e67896bedb756 |
| isotropic-elasticity--f7eb90376a.json | 2e775087ccf0d796f9a898478089bcfa585b52fbb49c34b683dac3e2182e73e8 |
| j2_props--2feae9f158.json | 60e0276a69451afabfbd7ea8d3853f067fe86c4adc830f8bcde3e20032e834ff |
| lemaitredamagenonlocal--e7f8623ec9.json | 2fe9e6fa95694fcd3b5cc470e2d991609b48c033e4f6ff4c2359da5469fc0882 |
| neohookean_umat--10759f1ffd.json | c36a4800f068315874042de6fc4fd2842b41a5cf934748f2617fee4f728146bd |
| umat--b8fa38353e.json | 2f6f520e7134db08c08204d443d19632a9b9c8b561dee9626bbae810d32ff1cf |
| umat_biofilm_visco_phase2--42eacbd5e8.json | dca7ea1336133056914c0ee9cd983c5bfe5e7aad523ef74761a308a2fa9045d6 |
| umat_elastic--7e9bb4c291.json | d5f4bb6c8ec86cc845678cba0f3fd53ba5593e36d1267f4ed5238f117b6bca2c |
| umat_iso_stretch--890619c18d.json | e79259989183d14e3977558fdb1d0cffc32da66465274b38964df1d621647a98 |
| umat_viscoelastic--c6ae96a734.json | 0e23f6c751fe90d5aab6a1c7ec97bde37c395c9aa392da86a11eeb074c840f41 |

## Commands And Checks

Environment used (confirmed imports printed both recovery paths and
/home/ammslab3/otilib/build_py311/pyoti/sparse.cpython-311-x86_64-linux-gnu.so):

```bash
WORKSPACE=/home/ammslab3/softwarex_work
RA="$WORKSPACE/imq-ra-recovery"
UMAT="$WORKSPACE/imq-umat-recovery"
PY="$WORKSPACE/.venv/bin/python"
cd "$RA"
export PYTHONPATH="$RA:$UMAT/src:/home/ammslab3/otilib/build_py311"
export UMAT_OTI_REPO="$UMAT"
export PYOTI_PATH=/home/ammslab3/otilib/build_py311
export OTILIB_ROOT=/home/ammslab3/otilib/build_py311
export RUN_OTILIB_TESTS=1
source /opt/intel/oneapi/compiler/latest/env/vars.sh
```

Actual regeneration commands, in order. The first ran before the archive move;
for reproduction use its now-archived input path. The tool refuses existing
--work directories, so reproduce into a new directory below recovery_fixtures.

```bash
"$PY" scripts/regenerate_recovery_fixture.py --umat "$UMAT" --cache "$WORKSPACE/discovery_cache" --results "$WORKSPACE/corpus_run/_superseded/pass12/results/store_verification.jsonl" --fixture tests/fixtures/verified/isotropic-elasticity--f7eb90376a.json --work "$WORKSPACE/imq_abaqus/recovery_fixtures/elasticity"
"$PY" scripts/regenerate_recovery_fixture.py --umat "$UMAT" --cache "$WORKSPACE/discovery_cache" --results "$WORKSPACE/corpus_run/_superseded/pass12/results/store_verification.jsonl" --fixture tests/fixtures/historical/94a92c01814f107a/j2_props--2feae9f158.json --work "$WORKSPACE/imq_abaqus/recovery_fixtures/j2"
"$PY" -m residual_core.core.fixture_residual_check --fixture tests/fixtures/verified/isotropic-elasticity--f7eb90376a.json --out docs/evidence/recovery_fixture_elasticity_check.json
"$PY" -m residual_core.core.fixture_residual_check --fixture tests/fixtures/verified/j2_props--2feae9f158.json --out docs/evidence/recovery_fixture_j2_check.json
"$PY" -m pytest -q tests/contract/test_a_fixture_is_refused_at_the_wrong_fingerprint.py --tb=short --junitxml=docs/evidence/recovery_fixtures_focused.xml
"$PY" -m pytest -q tests/contract tests/framework/test_a_fixture_is_held_to_the_rule_that_froze_it.py tests/framework/test_a_step_is_not_decoration.py tests/framework/test_a_verified_deck_drives_the_global_assembly.py tests/framework/test_a_verified_umat_drives_the_assembler.py tests/framework/test_the_assembled_residual_derivative_is_differenced.py tests/framework/test_the_assembled_state_derivative_is_differenced.py tests/framework/test_which_layer_the_failure_came_from.py -m 'not abaqus and not arc and not network' --tb=short --junitxml=docs/evidence/recovery_fixtures_subsystem.xml
"$PY" -m pytest -q -ra -m 'not abaqus and not arc and not network' --tb=short --junitxml=docs/evidence/recovery_fixtures_offline.xml
git diff --check
```

| Check | Result |
| --- | --- |
| First regression before regeneration | 1 failed (no current fixture), 1 passed (historical rejection) |
| First installed current fixture | 1 passed |
| Historical/current split, focused numerical checks | 59 passed |
| Final focused archive/current CLI checks | 21 passed, 2.89 s |
| Fixture/contract subsystem | 190 passed, 11 existing skips, 13.91 s |
| Full RA offline, no collection bypass | 419 passed, 11 existing skips, 57.05 s |

The 11 skips are unchanged pass11 cross-reader tests looking at the old
corpus_run/pass11 path. Targeted lookup found the superseded pass11 report under
corpus_run/_superseded/pass11; this task did not redirect those unrelated tests
or claim that their comparisons ran. No new skip was introduced. All 30 formerly
uncollected deck-module tests now execute, and four current/archive regressions
were added, explaining the increased denominator.

No blocker remains for this bounded suite repair. This is not a whole-corpus
verification, a current ten-material claim, a new wheel/clean-clone gate, or
proof of generic finite-strain/stateful UMAT support.