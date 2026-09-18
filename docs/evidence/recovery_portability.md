# Recovery Portability Slice

Date: 2026-09-18. **Bounded working-tree portability slice: PASS.** This is
not committed clean-clone acceptance, full historical corpus coverage, or
whole-ledger completion. The exact records, all original failure identities,
final skip identities, wheel hashes, import paths and raw-artifact digests are
in [recovery_portability.json](recovery_portability.json). The earlier
[clean-clone failure report](recovery_clean_clone.md) remains historical.

## Results

| Check | Result |
| --- | --- |
| Full RA offline suite | 437 passed, 11 skipped, 0 failures/errors |
| Full UMAT offline suite | 3249 passed, 171 skipped, 5 deselected, 0 failures/errors |
| Ten installed examples R-X1 through R-X5 and U-X1 through U-X5 | All passed, sequential execution |
| Fresh empty environment, newly built wheels, installation gate | Passed; no editable installs |
| Installed UMAT browser demo and same-contract CLI | Transform and compile passed; 21/21 structural checks, zero blockers |
| Browser desktop/mobile | 1440x1100 and 390x844; no page errors or horizontal overflow; server stopped |
| Existing genuine ODB, installed source-denied request | Passed independent analytic checks; no new Abaqus job |
| Producer fingerprint before/after | Both `6aa20d22e37f14c9` |

Raw work is retained at `/tmp/imq_portability_3W4QbqeM`. The source-suite
snapshots are `consumer-arbitrary` and `producer-arbitrary`, not canonical or
recovery folder names. They were local clones with tracked/untracked working
changes overlaid, not clean committed inputs. The successful installation
gate is under `installed-final`; its conservative `wheel-from-working-tree`
label is retained. Test-only locator changes were then synced to the producer
snapshot and its whole suite rerun. The later presentation consume-script
repair was checked directly against the installed gate artifacts.

| Wheel | SHA-256 |
| --- | --- |
| residual_assembler-0.1.0-py3-none-any.whl | `2dd7ece45b91d14249615cdbb184d627b914a10809d41179894161c5ba24d5ff` |
| umat_oti-1.1.0-py3-none-any.whl | `33a7f6ae8819b97cef88af96bcab9cebdb57883206158aca3df7cce25755f98b` |

## Failure Classification

The JSON classifies **every one of the 114 failure/error entries** in the two
original attempts, retaining run and test identity. Repeated failures are not
counted as distinct defects:

| Original cause | Entries across both attempts | Resolution |
| --- | ---: | --- |
| RA hardcoded producer sibling | 98 setup errors + 2 failures | Explicit `UMAT_OTI_REPO`, canonical default, shared test helper; reproducer accepts any cwd |
| UMAT browser prerequisite | 4 errors + 2 failures | Explicit compatible Playwright 1.48.0 and Chromium 1140 |
| UMAT compiler prerequisite | 2 failures | Explicit ifort PATH and compiler runtime libraries |
| UMAT repository standards | 2 failures | Portable live command; exact historical transcript exemptions with reasons |
| UMAT documentation links | 2 failures | Canonical companion links, explicit `RESASM_REPO` resolution, missing cache artifacts described truthfully |
| Historical source discovery | 2 failures | Existing skip now checks the actual required original files, not merely their parent directory |

Other reported workflow defects were also handled: the existing uncommitted
RA template fix was preserved (all original 11 tests plus two locator tests
passed); the UMAT wheel now carries only the original public demo contract
and its source; R-X4 defaults to isolated installed backend imports; and the
example auditor honors supplied OTILib selectors rather than developer HOME.
The presentation consumer no longer injects project sources. Its optional
licensed preparation retains its existing job-directory safety restrictions
and was not executed.

## Implementation Boundaries

UMAT packaging stages two byte-identical originals into app resources while
preserving their relative source reference. `importlib.resources` resolves
them into a writable workspace. GUI and CLI still use the same transformation
service; no second mathematical implementation or whole-repository copy was
introduced. Wheel resource, arbitrary-path and real compile regressions pass.

The existing fingerprint rule already excludes app code; no exemption was
added or widened. No fingerprinted transform code, generation/lock file or
fixture bytes changed. There is therefore no new generation to retire or
regenerate. Test extras now declare wheel, mpmath and matplotlib. SciPy is an
explicit external OTILib runtime prerequisite, not a replacement for OTILib.

Cross-project fixture tests use explicit checkout helpers. Broken explicit
configuration raises instead of silently selecting another checkout. The
documentation auditor still rejects missing companion targets; its regression
tests cover arbitrary names and broken links. Historical transcript exemptions
are exact named files, never source/tests or all evidence/documents.

## Reproduction Commands

Set `RA`, `UMAT`, `BASE_PYTHON`, `ODB`, `OTILIB_BUILD`, `INTEL_BIN`,
`INTEL_LIBS` and `BROWSERS` to explicit local inputs. `ODB` must be the existing
genuine presentation J2 input, not synthetic data. The run used OTILib built
at pinned upstream `a4b7a05ca275e8d441b0b717b7b272e96728ffcc` in the earlier
clean-clone investigation, not a hidden developer build. Its SciPy prerequisite
was installed into the new wheel environment. Browser and compiler versions
are external prerequisites, not project-source import paths.

```sh
WORK=$(mktemp -d /tmp/imq_portability_XXXXXXXX)
"$BASE_PYTHON" "$RA/scripts/clean_install_gate.py" \
  --ra-repo "$RA" --umat-repo "$UMAT" --python "$BASE_PYTHON" \
  --odb "$ODB" --work "$WORK/installed-final"
PY="$WORK/installed-final/env/bin/python"
"$PY" -I -m pip install 'scipy>=1.10' 'playwright==1.48.0'
export HOME="$WORK/installed-final/home"
export UMAT_OTI_REPO="$UMAT" RESASM_REPO="$RA"
export PYOTI_PATH="$OTILIB_BUILD" OTILIB_ROOT="$OTILIB_BUILD"
export RUN_OTILIB_TESTS=1 PYTHONNOUSERSITE=1
export PLAYWRIGHT_BROWSERS_PATH="$BROWSERS"
export PATH="$(dirname "$PY"):$INTEL_BIN:/usr/local/bin:/usr/bin:/bin"
export LD_LIBRARY_PATH="$INTEL_LIBS"
cd "$RA"
PYTHONPATH="$RA:$UMAT/src:$OTILIB_BUILD" "$PY" -m pytest -q -ra \
  -m 'not abaqus and not arc and not network' --junitxml="$WORK/ra_offline_final.xml"
cd "$UMAT"
UMAT_OTI_SNAPSHOT_ROOT="$RA/sources" PYTHONPATH="$UMAT/src:$OTILIB_BUILD" \
  "$PY" -m pytest -q -ra -m 'not abaqus and not arc and not network and not corpus_pass' \
  --junitxml="$WORK/umat_offline_final.xml"
unset PYTHONPATH PYTHONHOME
"$PY" -I "$RA/scripts/audit_recovery_usage.py" --umat "$UMAT" \
  --phase examples --imports installed --work "$WORK/examples" \
  --evidence-dir "$WORK/example-evidence"
"$PY" -I "$RA/scripts/reproduce_presentation_request.py" consume \
  --work "$WORK/installed-final" --out portability_consumed
"$PY" -I -m pip check
```

Actual isolated runs used `env -i`, scratch HOME and explicit PATH. Source
suites intentionally imported the source snapshots; installed examples used
`-I` and no project PYTHONPATH. R-X4 records both backend paths in site-packages
and reads model inputs explicitly from the snapshots. Browser harness source
is embedded in JSON; it unsets both UMAT checkout selectors, launches the
installed app, clicks the real demo button and invokes the installed CLI on
the same extracted contract. All owned servers are terminated and waited.

## Remaining Evidence And Handoff

RA's 11 skips are the absent historical pass11 frozen store. UMAT's 171 skips
name absent pass9/pass10/pass11 runs, external source/cache/deck inputs and
other existing unavailable historical cases. All identities/reasons are in
JSON. Five licensed corpus tests remain deselected by the documented offline
selection. No large external corpus was copied, no missing field assertion
was relaxed, and no synthetic integration substituted for missing originals.
The available documented permissive submodule was included explicitly.

Intermediate issues are retained: a copied submodule Git pointer was removed
only from the temporary snapshot; the first fresh RA environment lacked
SciPy; and one browser harness CLI call omitted `--config` before correction.
The second UMAT suite was requested synchronously but auto-backgrounded by
the terminal tool after output idle; later browser/presentation checks
overlapped it. The ten examples themselves ran sequentially. Final suite
results were collected only after completion.

No commits, pushes, recovery-branch changes, subagents or new Abaqus jobs were
performed. **Next lead action:** review and commit the preserved template fix,
these portability changes and evidence in both recovery repositories; then
run a new exact-commit clean-clone gate with explicitly supplied external
inputs. This report does not promote the entire ledger or claim missing
historical corpus evidence is complete.