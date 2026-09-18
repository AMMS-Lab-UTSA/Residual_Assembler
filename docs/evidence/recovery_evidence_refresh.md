# Recovery Evidence Refresh

Later fixture recovery: see [recovery_fixtures.md](recovery_fixtures.md) for two
genuine current reruns and the exact relocation map of all ten original RA
fixtures. This report and recovery_evidence_inventory.json remain the original
retirement snapshot, not an inventory of the later operational directory.

Date: 2026-09-18. The two UMAT offline failures were repaired with necessary
shared-contract synchronization into this RA worktree. No commits, pushes,
branch changes, subagents, Abaqus solves or corpus jobs. All prior work remains.
This is NOT a committed clean-clone gate, a wheel rebuild, or whole-ledger PASS.
Earlier wheel hashes in recovery_install.md still identify earlier builds.

## Exact Before And After

Both original UMAT tests reproduced before editing: **2 failed in 0.28s**.
`test_the_recorded_generation_is_this_worktrees_actual_transform` found recorded
94a92c01814f107a versus actual 6aa20d22e37f14c9.
`test_every_terminal_state_in_the_report_is_marked_external_or_internal` found
arguments_diverged_before_the_routine incorrectly EXTERNAL.

History inspected: d3bffa2 genuinely reran pass12 before recording its generation;
bcc8605 made divergence INTERNAL because the solver's later arguments depend on
each build's earlier outputs. The retained report had not been regenerated.
Shared brief, implementation audit and both installation records were read.

The current generation now records the measured transform 6aa20d22e37f14c9,
explicitly retiring older evidence. No old measurement received a new generation.
Contract version remains 3.0.0. UMAT's existing write_lock() produced combined
digest e326b594c40f95ebd1aa3ac7455d2060035b4d5faac45fafa97591bcc239a10e.
Only generation and lock were copied into schemas, byte-for-byte, as required
by transform_generation.json's how_to_update instructions. No ABI asset changed.

The canonical UMAT registry builder uses retained Record data, kind_of,
summarise and markdown through --refresh-retained. Exactly three row kinds move
external to internal. All other row fields, original generated date, numerical
observations, evidence/source references and historical fingerprint remain.
The separate evidence_currency metadata and report notice retire current claims.
verification_is_current inside retained rows means current relative to the OLD
store, not current code.

| Historical registry quantity | Before | After |
| --- | ---: | ---: |
| D1 acquired | 391 | 391 |
| D2 adequately specified | 260 | 260 |
| Fully verified / all gates true | 44 / 44 | 44 / 44 |
| External / internal | 122 / 225 | 119 / 228 |
| Divergence records (2 in D2) | 3 external | 3 internal |

## Exhaustive Retirement Inventory

[The shared inventory](recovery_evidence_inventory.json) enumerates 757 UMAT
files and 10 RA files by exact path, SHA-256, byte length, embedded generations
and parse error. Every supporting file in the listed historical evidence trees
is included, even without a generation field; hashes are not verification.

| Scope | Files | Retired current claims |
| --- | ---: | --- |
| UMAT umat tree | 313 | All 44 contracts at 94a92c01814f107a, baseline, registry, withdrawals, source identities, decks, histories and results |
| UMAT tests/fixtures tree | 23 | Four verified fixtures and both corpus bundles, including logs/decks/histories; corpus_current is a historical directory name |
| UMAT paper_results tree | 421 | All retained corpus reports, figures, tables, presentations and parameter-sensitivity results |
| RA tests/fixtures/verified tree | 10 | All pass12 fixtures retain original generation and original/converted numbers |

Older generations ff94800b1884bcc0, b0d27ee53c630500, b5c7a8d71b01c6e3 and
e4257779bd847cc4 are retained too. The refreshed registry's current fingerprint
is retirement metadata, NOT the generation of its observations. Two archived
ARC stdout files named .json are not JSON; hashes and parse errors are recorded.
No historical 57-terminal/44-gate pass12 claim, 18/18 benchmark claim, 20-model
sensitivity claim or old presentation proves current corpus capability. External
corpus_run/transform_store data are untouched and remain historical by generation.

Consumers inspected: UMAT contract.require_current/check_fingerprint and
TransformStore.current_entries/stale_entries; RA contract_reader.require_current,
materials.verified_fixture.load/load_all and core.fixture_residual_check.
The RA production loader and independent contract reader both reject all ten
real fixtures by name with both fingerprints. Original/converted arrays remain
readable. STORE_PROVENANCE keeps pass12's 244/57/44 counts with ORIGINAL
94a92c01814f107a, not the new dynamic generation. No bypass or skip was added.
Legacy framework helpers demanding current fixtures will refuse historical ones;
historical explicit-generation comparisons are not current capability evidence.

## New Bounded J2 Run

The current UMAT provider J2 test passed (1 test, 7.23s), followed by a fresh
verification CLI run retained in UMAT's
.pytest_cache/recovery_evidence_j2_run/verification.json. Current transformed
code and a separately compiled ORIGINAL UMAT were compared over seven increments:
elastic, elastic, plastic, plastic, elastic, plastic, elastic.
Comparison counts: 49 eval primal, 6 final march primal, 588 eval/504 march
parameter FD, 756 eval/756 march tangent FD. Worst parameter error
2.493046730365579e-7; worst tangent error 2.838732044089759e-8; unchanged 2e-6
tolerance. Primal stress/state errors: 5.684341886080802e-14 / 2.168404344971009e-19.
FD plateau stress/state/tangent: 1.994408000705079e-8 / 7.515721719070879e-10 /
2.7870660023874207e-10. Carry reset changes derivatives by 2720.125996569597.
Source hash: 9b779f0c6cadf9c4. Fresh object SHA-256:
8cb54ae11c1ecf511c4caaad7e05e7dc38a0e709ef9b86246d73bfcbb983d28d.
This is a real compiler/FD material-point run, NOT an Abaqus re-freeze of the
old four-step fixture. clean_install_verified=false. Generic UMAT/FCC,
full-sized models, publication/corpus reproduction and clean-clone packaging
remain unestablished by it.

## Measured Checks

| Check | Result |
| --- | --- |
| Original two UMAT failures before edits | 2 failed, 0.28s |
| New RA retirement regression before metadata repair | 1 failed, old fixture accepted |
| UMAT generation/lock after repair | 33 passed |
| Report module plus original generation guard | 19 passed |
| UMAT final focused contracts/report/documentation | 62 passed, 4.25s |
| RA actual-fixture/production-refusal module | 17 passed |
| RA contract/integration with OTILib path omitted | 148 passed, 2 failed, 11 skipped |
| RA contract/integration with documented OTILib environment | 150 passed, 11 skipped, 30.75s |
| Final RA contract/integration after all production edits | 150 passed, 11 skipped, 31.23s |
| Full final UMAT offline suite | 3324 passed, 0 failed, 125 skipped, 5 deselected, 6 warnings, 477.24s |

Final JUnit: UMAT .pytest_cache/recovery_evidence_offline.xml (3449 executed,
zero failures/errors); RA .pytest_cache/recovery_evidence_contract_integration_final.xml.
UMAT warnings: four existing deprecated transform-entrypoint calls and two
constrained-layout warnings in a figure-rejection test. Its 125 existing skips
remain unavailable historical evidence/previously documented gates. Actual code
still matches the recorded generation. Both shared asset copies match. A Git
HEAD comparison proves exactly three record kind fields changed; no historical
fixture or material numerical evidence file changed, only the collection notice.

RA's 11 skips are existing unavailable pass11 cross-reader records. The failed
environment attempt is kept in separate JUnit. No tests deleted, tolerances
weakened, hidden skips or new xfails. RA full legacy framework suite is not
represented by the contract/integration count.
Current-generation legacy-fixture workflows require a real new Abaqus fixture
run; the old fixtures are correctly refused. No full legacy-framework PASS,
clean-clone gate, regenerated provider wheel hashes or whole-ledger completion
is claimed.

## Reproduce

```sh
WORKSPACE="$HOME/softwarex_work"
UMAT="$WORKSPACE/imq-umat-recovery"
RA="$WORKSPACE/imq-ra-recovery"
PY="$WORKSPACE/.venv/bin/python"
cd "$UMAT"
export PYTHONPATH="$UMAT/src"
"$PY" -c 'from umat_oti.store import transform_fingerprint; print(transform_fingerprint())'
"$PY" -c 'from umat_oti.contract import write_lock, verify_lock; print(write_lock()); verify_lock()'
cp src/umat_oti/contract/schemas/transform_generation.json "$RA/schemas/transform_generation.json"
cp src/umat_oti/contract/contract_lock.json "$RA/schemas/contract_lock.json"
"$PY" tools/build_corpus_registry.py --refresh-retained paper_results/corpus/corpus_registry.json
"$PY" tools/inventory_retired_evidence.py --umat "$UMAT" --ra "$RA" --out docs/evidence/recovery_evidence_inventory.json
cp docs/evidence/recovery_evidence_inventory.json "$RA/docs/evidence/recovery_evidence_inventory.json"
"$PY" -m pytest -q tests/test_contract_fixtures.py tests/test_contract_schema_lock.py tests/test_the_two_denominators_stay_apart.py tests/test_repository_standards.py
"$PY" -m pytest -q tests/test_provider_recovery.py::test_provider_j2_original_fd_and_carryover --junitxml=.pytest_cache/recovery_evidence_j2.xml
"$PY" -m umat_oti.validation.parameter_sensitivity_provider parameter_sensitivity/models/m3_j2/contract_v2.json --out .pytest_cache/recovery_evidence_j2_run
"$PY" -m pytest -q -ra -m 'not abaqus and not arc and not network and not corpus_pass' --junitxml=.pytest_cache/recovery_evidence_offline.xml
cd "$RA"
env PYTHONPATH="$RA:$UMAT/src:$HOME/otilib/build_py311" UMAT_OTI_REPO="$UMAT" PYOTI_PATH="$HOME/otilib/build_py311" OTILIB_ROOT="$HOME/otilib/build_py311" RUN_OTILIB_TESTS=1 "$PY" -m pytest -q -ra tests/contract tests/integration -m 'not abaqus and not arc and not network' --junitxml=.pytest_cache/recovery_evidence_contract_integration_otilib.xml
```

All runs use explicit interpreter, cwd and recovery import roots and run
synchronously/sequentially. No detached one-shot process or inherited main-tree
imports. Future source-generation changes require another honest retirement or
executable re-freeze, never changing old fixture IDs to imply currency.