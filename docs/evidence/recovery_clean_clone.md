# Committed Clean-Clone Verification

Date: 2026-09-18. **Full requested gate: FAIL / incomplete.** Bounded wheel
installation and the genuine J2 presentation workflow passed. Neither the
whole offline suites nor every installed example/GUI workflow passed.
The complete structured record is [recovery_clean_clone.json](recovery_clean_clone.json).

## Exact Inputs And Isolation

Both LOCAL recovery repositories were cloned with `--no-hardlinks --branch
integration/imqcam-recovery-2026-09-18` into the new parent
`/tmp/imqc_clean_clone_20260918_RwJ8VT`:

| Canonical clone | Full committed SHA |
| --- | --- |
| UMAT_source_transformation | dbd8ff229023b57e4231a2250a54e97dc111f547 |
| Residual_Assembler | 5ea3da4bdd250fe6eea79677cec44e59e86de197 |

Source checkouts and both clones had empty porcelain status before validation.
Both full-suite attempts also recorded empty pre-test status. No clone source
or test was edited. The final unmodified example auditor generated changes to
only `docs/evidence/usage_examples.json` in each clone, after the suites; those
are not committed inputs. No commits, pushes, new branches, or subagents.

The unchanged UMAT forwarding gate was run from its clone, with explicit RA
clone, `--python /home/ammslab3/anaconda3/bin/python3.11`, the documented
external ODB below, and new `--work /tmp/imqc_clean_clone_20260918_RwJ8VT/gate`.
It created an empty Python 3.11.7 venv, built both wheels, and installed RA's
`gui,test,yaml` and UMAT's `test` extras. All 20 child commands exited zero.
The gate's conservative raw `wheel-from-working-tree` label is preserved in
JSON; the separately recorded clean clone provenance establishes this exact
committed pair, not a later follow-up revision.

Installed commands used `env -i`, scratch HOME, no editable installs, no
inherited PYTHONPATH, and no user site. `-I` probes resolved project modules
only under the new venv's site-packages. Source suites deliberately used the
cloned source trees and a freshly built external OTILib; these are separate
from wheel-only execution. Base standard-library paths are not source leaks.
Final `pip check` passed. Full dependency versions and every recorded command,
cwd, exit and raw-log digest are in JSON. Raw logs, temporary harness sources,
JUnit and screenshots remain under the temporary parent; harness sources are
also embedded in JSON. Early read-only routing failures are disclosed there.

| Milestone wheel | SHA-256 |
| --- | --- |
| residual_assembler-0.1.0-py3-none-any.whl | 8c38bb185fd162e336e20ca6711f54bafd3d004e13e477e9c452bcbcc39bf4e2 |
| umat_oti-1.1.0-py3-none-any.whl | f0cf7f51b59493c7d8cbd0cd7e8030b3ead56c4071e5ad192d077d87b8ef053d |

## Fresh External Dependencies

OTILib was cloned afresh from `https://github.com/mauriaristi/otilib.git`,
commit `a4b7a05ca275e8d441b0b717b7b272e96728ffcc`, GPLv3. No existing external
build or direction table was used. The committed setup helper exited 1 because
it requires Conda. A no-Conda build succeeded using CMake 3.31.10, Cython 3.3.0,
NumPy 2.4.6, SciPy 1.17.1 and GCC/gfortran 9.4.0:

```sh
cmake -S "$CLEAN_ROOT/otilib" -B "$CLEAN_ROOT/otilib-build" -DBUILD_TESTING=OFF
cmake --build "$CLEAN_ROOT/otilib-build" --target oticython --parallel 2
cmake --build "$CLEAN_ROOT/otilib-build" --target gendata --parallel 2
```

All three exited zero; genuine sparse first/second derivative assertions
passed. Both explicit OTILib selectors pointed to this new build directory.
Upstream OTILib's own test suite was not run. The proven procedure is now a
separate follow-up in [OTILIB_VENV.md](../OTILIB_VENV.md), not in the milestone.
RA's documented permissive-source setup fetched only the pinned Oxford
Crystal Plasticity submodule; restricted tiers stayed uninitialized.

All executions were requested synchronously and sequentially. The terminal
tool automatically moved the long OTILib command to a background handle after
output idle; it completed exit zero before the next build/test execution.
There was no intentional detached command. Browser servers were owned by
synchronous harnesses and terminated/waited in `finally` blocks.

## Real ODB And Presentation

The explicitly declared existing scientific input was
`imq_abaqus/recovery_presentation/imqrp_reference/collaborator/Analysis.odb`,
SHA-256 `54745b97b515de159c3251e8a3ba0ed89ea1c65452e388ef34c8c81794e3b736`.
It was not treated as an undeclared development cache.

One new sequential Abaqus job, `imqc_j2`, completed under
`imq_abaqus/recovery/imqc_clean_clone_RwJ8VT`. Only the public deck and J2
material source were copied from the clones. Toolchain: Abaqus 2021.HF5,
ifort 2021.10.0; Standard checked out five licensed tokens. The committed
prepare helper hardcodes `imqrp_` and `imq-umat-recovery`, so its documented
public job recipe was executed directly with the required `imqc_` prefix;
no aliases or patched clone scripts hid that limitation.

Fresh ODB SHA-256:
`92460c0ec883a14e0ea4a96021e52cb471dd509f93a4c5dbf0918ea9c9a3c54e`.
Public deck/source hashes and license-file hashes are in JSON. Both projects
declare GPL-3.0-only; OTILib is GPLv3; Abaqus is proprietary licensed software.
No independent ODB redistribution grant was established and no ODB/binary was
added to either repository.

The installed milestone `resasm request` consumed both ODBs successfully.
Four increments and four scalar results were checked. Fresh-ODB numerical
results exactly matched the existing-ODB results. Deleting ONLY the gate's
new `fresh_odb_results` and regenerating produced byte-identical JSON, CSV and
text reports. The original gate's source-denied replay also matched all three
public outputs. Ordinary output correctly retains `verified=false`.

Independent uniaxial reference: `U1=300/E+(300-SIGY0)/H`. Relative derivative
errors were E `1.6568627353483993e-7`, SIGY0 `8.673617379884035e-16`, and H
`8.839354487082118e-9`, all below `2e-5`; nu absolute error was below `1e-8`.
Reaction and stress errors were below `1e-3`, state error below `1e-6`.

## Examples And Browsers

All ten auditor cases were attempted twice. Only this gate's sample output
directory was inventoried, deleted and regenerated between rounds. Eight
passed twice with identical scientific reports: R-X3, R-X4, R-X5 and U-X1
through U-X5. R-X4's committed reproducer injects clone sources despite the
explicit canonical `--provider-repo`; it is NOT installed-command proof.
Thus only seven cases have successful installed-module evidence at these SHAs.
R-X1/R-X2 failed in both rounds: the milestone RA wheel omits `templates/`.
The unmodified auditor itself exited 1 before examples because it overwrites
the supplied OTILib selectors with `$HOME/otilib/build_py311`. No historical
scratch directory was substituted to make it pass.

Playwright 1.63.0 could not download Chromium for Ubuntu 20.04. Version 1.48.0
then downloaded fresh Chromium 130.0.6723.31 (build 1140) into scratch HOME.
Real servers launched both apps by installed module path. RA's actual browser
request passed, all three downloads matched generated files, and numerical
results matched CLI results. Both apps rendered at 1440x1100 and 390x844 with
no exception or horizontal overflow; all servers stopped. UMAT displayed
"The demo contract is missing" at an erroneous environment-relative path.
Its installed transformation workflow therefore remains unverified, despite
successful startup. No all-GUI-workflows PASS is claimed.

## Full Offline Suites

| Run | Passed | Failed | Errors | Skipped | Deselected |
| --- | ---: | ---: | ---: | ---: | ---: |
| RA initial | 372 | 1 | 49 | 19 | 0 |
| RA with documented permissive source | 380 | 1 | 49 | 11 | 0 |
| UMAT initial | 3234 | 7 | 4 | 172 | 5 |
| UMAT with explicit ifort and compatible browser | 3243 | 3 | 0 | 171 | 5 |

RA used `-m 'not abaqus and not arc and not network'`; UMAT additionally
excluded `corpus_pass`, as documented. Neither suite selected only passing
tests. Final runtimes: RA 28.203 s, UMAT 582.051 s. These are not full online,
licensed corpus, or all-platform suites; existing skips remain incomplete.

RA's 49 errors and one failure all reference nonexistent sibling
`imq-umat-recovery`, despite canonical folder names and UMAT_OTI_REPO. UMAT's
three remaining failures are repository standards (absolute development
paths), documentation links (recovery-specific sibling paths), and source
discovery (no corpus source available). No symlink, copied developer cache,
new skip, weakened tolerance or relabelled fixture was used to hide them.

## Separate Follow-Up

A concrete RA packaging fix was made ONLY in the original recovery checkout:
ship the five existing public templates as wheel data and discover their
installed location. A separate follow-up wheel/environment successfully
created and ran the two previously failing templates; all five analytic
derivative directions in each matched within `1e-8`. Focused regressions:
**11 passed, zero failures/skips**. Follow-up wheel SHA-256:
`7dbb43f95426529742f46fff09fffcf505e635c615946583804b4c61895c87f0`.
This is NOT proof for RA 5ea3da4 and requires a new commit and a new clean clone.
The milestone environment was never replaced by the follow-up wheel.

Only U-A1, U-A4, R-A1 and R-A6 are promoted in each completion ledger, bounded
to exact-SHA wheel installation and installed CLI startup. All other rows stay
unchanged; CI/full-suite, every-example, every-GUI, no-path-dependency and
whole-ledger completion are explicitly NOT established. The historical
usage index remains a previous audit snapshot; these four row-specific proofs
are in this report and JSON. Evidence and fixes are uncommitted follow-up work.