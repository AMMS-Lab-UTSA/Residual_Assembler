# Recovery Wheel Installation Evidence

**Historical installation snapshot:** the two UMAT failures recorded below
have a subsequent [evidence repair record](recovery_evidence_refresh.md).
These wheel hashes still identify the earlier builds; they were not refreshed
or promoted to a committed clean-clone gate by that repair.

Date: 2026-09-18. Both recovery worktrees remained on
`integration/recovery-2026-09-18`. Existing changes were preserved.
No commits, staging, pushes, branch changes, old-worktree edits, VS Code config
edits, subagents, detached one-shot tests, or new Abaqus analysis jobs.

## Result And Boundary

**PASS: wheels built from the current working trees, installed into a new venv
outside both repositories with a scratch HOME.** This is not a committed
final-branch clean-clone gate. The repositories are still dirty and contain
untracked implementation files. The completion ledger is not declared complete.

The pre-existing packaging changes (provider console registration and RA replay
contract package data) work as written. No additional runtime dependency or
plotting extra was needed for this bounded workflow. RA uses its existing
`gui,test,yaml` extras; UMAT uses `test` plus its declared runtime dependencies.
Publication plotting remains in UMAT's existing `paper` extra. This gate does
not select RA's `bridge` extra, which points at a different historical revision;
it installs the companion wheel built alongside RA instead.

## Supported Installation

Use healthy Python >=3.10 for the joint workflow, with `venv`, `ssl`, and
`ctypes`. Linux/gfortran is the tested platform. From an external directory,
with both source checkouts available for building:

```sh
python3.11 -I -m venv /tmp/imqr_user_env
/tmp/imqr_user_env/bin/python -I -m pip wheel --no-deps --wheel-dir /tmp/imqr_wheels /path/to/imq-ra-recovery
/tmp/imqr_user_env/bin/python -I -m pip wheel --no-deps --wheel-dir /tmp/imqr_wheels /path/to/imq-umat-recovery
/tmp/imqr_user_env/bin/python -I -m pip install '/tmp/imqr_wheels/residual_assembler-0.1.0-py3-none-any.whl[gui,test,yaml]' '/tmp/imqr_wheels/umat_oti-1.1.0-py3-none-any.whl[test]'
/tmp/imqr_user_env/bin/python -I -m pip check
```

Use new directories. Network access is needed to resolve/download build and
runtime dependencies unless already supplied through an appropriate wheelhouse.
No editable install, PYTHONPATH or checkout working directory is needed to run:

```sh
/tmp/imqr_user_env/bin/umat-oti-provider build /path/to/public_model/contract_v2.json --out /tmp/new_provider
cd /path/to/shared_artifacts
/tmp/imqr_user_env/bin/resasm request --model Analysis.inp --odb Analysis.odb --material OTI_UMAT.obj --request sensitivity_request.json --out new_results
```

Only the developer runs the first command with public source/contract inputs.
The collaborator directory contains `Analysis.inp`, `Analysis.odb`,
`OTI_UMAT.obj`, `Mapping.json`, and `sensitivity_request.json`. Public outputs
are exactly `sensitivity_results.json`, `sensitivity_tables.csv`, `run_report.txt`;
full fields and link products stay under `private/`. Ordinary results retain
`verified=false`; the gate's independent analytic check is separate.

Installed GUI launch, with no checkout launcher or sys.path injection:

```sh
PY=/tmp/imqr_user_env/bin/python
RA_APP=$("$PY" -I -c 'import importlib.util; print(importlib.util.find_spec("residual_core.app.streamlit_app").origin)')
UMAT_APP=$("$PY" -I -c 'import importlib.util; print(importlib.util.find_spec("umat_oti.app.streamlit_app").origin)')
"$PY" -I -m streamlit run "$RA_APP" --server.address=127.0.0.1 --server.port=8501
# In a separate terminal:
"$PY" -I -m streamlit run "$UMAT_APP" --server.address=127.0.0.1 --server.port=8502
```

Choose free ports. UMAT's existing GUI is not a newly implemented provider-build
GUI. Repository-only example discovery and publication/corpus workflows are not
established by an installed launch check.

## Exact Executed Gate

The final run was invoked through the UMAT forwarding script:

```sh
cd /home/ammslab3/softwarex_work/imq-umat-recovery
/home/ammslab3/softwarex_work/.venv/bin/python scripts/clean_install_gate.py \
  --ra-repo /home/ammslab3/softwarex_work/imq-ra-recovery \
  --python /home/ammslab3/anaconda3/bin/python3.11 \
  --odb /home/ammslab3/softwarex_work/imq_abaqus/recovery_presentation/imqrp_reference/collaborator/Analysis.odb \
  --work /tmp/imqr_install_20260918_03
```

For a repeat, choose a NEW `--work` path or omit it for `mkdtemp`. The RA entry
is `scripts/clean_install_gate.py --umat-repo /path/to/imq-umat-recovery` with
the same `--python`, `--odb`, and optional `--work`. Explicit branch checks reject
the wrong worktrees. Missing executables, incompatible Python, missing ODB,
build/installation errors, failed assertions or calculation failures exit
nonzero; no core mocks, manufactured ODB, fallback calculation or new skips.

The shared venv is used only as bootstrap. Its base interpreter was discovered
and checked: `/home/ammslab3/anaconda3/bin/python3.11`, Python 3.11.7.
`/usr/local/bin/python3.11` was actually tested and fails importing `_ctypes`;
`/usr/bin/python3` is healthy but 3.8.10, below the joint package requirement.

The gate strips PYTHONPATH/PYTHONHOME, pip overrides, VIRTUAL_ENV, and OTILib
development selectors; disables user-site and pip config; uses
`/tmp/imqr_install_20260918_03/home` and its XDG config/cache. The installed
interpreter is `/tmp/imqr_install_20260918_03/env/bin/python`. The actual commands,
cwd, return codes, wheel/input hashes, resources, sys.path and all versions are
retained in `/tmp/imqr_install_20260918_03/report.json`, with logs under `logs/`.
Wheel building used pip's isolated build environment, not bootstrap packages.

Actual installation command after both `pip wheel --no-deps --wheel-dir ...`
commands and `pip install --upgrade pip` in the new environment:

```sh
/tmp/imqr_install_20260918_03/env/bin/python -I -m pip install '/tmp/imqr_install_20260918_03/wheels/residual_assembler-0.1.0-py3-none-any.whl[gui,test,yaml]' '/tmp/imqr_install_20260918_03/wheels/umat_oti-1.1.0-py3-none-any.whl[test]'
/tmp/imqr_install_20260918_03/env/bin/python -I -m pip check
```

`pip check`: exit 0, no broken requirements. Installed versions: Python 3.11.7,
pip 26.2.1, residual-assembler 0.1.0, umat-oti 1.1.0, numpy 2.4.6,
pandas 3.0.6, streamlit 1.64.0, sympy 1.14.0, pytest 9.1.1, PyYAML 6.0.3.
The full transitive version list is in the JSON report; dependency resolution
is not locked, so future resolution need not be byte-identical.

Wheel SHA-256:

| Wheel | SHA-256 |
| --- | --- |
| residual_assembler-0.1.0-py3-none-any.whl | 4eca0c2c2201fb707b7423c124344f40b29c84e57c75c9b3bd0c0c3f6770a5fd |
| umat_oti-1.1.0-py3-none-any.whl | ca067f9e762359639d573efd93f533e0095c167eff4f53ab3ac8691d55ca6881 |

## Installed Evidence

Every checked import is below
`/tmp/imqr_install_20260918_03/env/lib/python3.11/site-packages/`:
`residual_core/__init__.py`, `resasm_user/__init__.py`, `umat_oti/__init__.py`,
`umat_oti/provider/build.py`, `residual_core/replay/presentation.py`, and both
`residual_core/app/streamlit_app.py` and `umat_oti/app/streamlit_app.py`.
sys.path contains only base standard-library locations and this venv's
site-packages; neither recovery tree is present. No distribution is editable.

Resource assertions cover RA's example manifest, Fortran replay driver/stubs,
replay `CONTRACT_VERSION.json` and ABI header; UMAT's support Fortran sources,
OTI templates, contract schemas and lock. Fresh provider compilation also
exercises runtime support generation from installed resources.

`umat-oti-provider build` compiled copied public m3_j2 source/contract outside
the source tree. `resasm request` ran from the artifact-only collaborator
directory, extracted all genuine ODB frames through installed exporter code,
replayed all four increments and returned four scalar outputs. A second run of
the installed console script under a Python audit hook denied material Fortran
source reads except its exact generated link shim. All three public files were
byte-identical to the normal console run. This is an audited check, not an OS
sandbox; native compiler/Abaqus processes are not sandboxed.

Independent uniaxial J2 reference: `U1=300/E+(300-SIGY0)/H`.
Relative derivative errors: E `1.6568627353483993e-7`, SIGY0
`8.673617379884035e-16`, H `8.839354487082118e-9`; threshold `2e-5`.
The zero-reference nu derivative is checked with absolute threshold `1e-8`.
No tolerance was weakened. The genuine ODB is the previous licensed reference,
SHA-256 `54745b97b515de159c3251e8a3ba0ed89ea1c65452e388ef34c8c81794e3b736`.
Neither repository tracks an ODB; that external fixture remains required.

Both installed GUI modules rendered using `streamlit.testing.v1.AppTest` with
zero exceptions and nonempty titles/tabs. Separate real Streamlit servers
reached HTTP health readiness and were stopped/waited by exact PID:
RA PID 1065394, port 44163; UMAT PID 1065413, port 33311. Neither URL is live.
This run checks installed rendering and server launch, not browser downloads;
the prior source-tree browser proof remains in `recovery_presentation.md`.

Run history: `_01` failed an incorrect namespace-resource probe before provider
execution; the resource existed in the wheel, but converting a namespace
resource handle to a filesystem string was wrong. It was corrected to exact
paths under the owning package and immediately rechecked against that wheel.
`_02` passed the complete initial gate. `_03` passed the strengthened gate above.
Failed reports remain retained; no failure was relabelled as a pass.

## Regression Results

| Check | Passed | Failed | Skipped | Deselected |
| --- | ---: | ---: | ---: | ---: |
| New gate regressions across both repositories | 8 | 0 | 0 | 0 |
| Full relevant RA offline suite | 390 | 0 | 11 | 0 |
| Full relevant UMAT offline suite | 3316 | 2 | 125 | 5 |

RA runtime 53.39 s; UMAT 495.26 s. JUnit files are
`.pytest_cache/recovery_install_offline.xml` in each repository. Existing skips
require unavailable frozen corpus artifacts or other documented prerequisites;
none were introduced to make this work pass. Source-tree suites use the shared
development venv and are separate from the wheel-isolation evidence.

After the guide edits, UMAT `tests/test_repository_standards.py` plus
`tests/test_clean_install_gate.py` passed 10/10. Its existing
`tools/audit_documentation_commands.py --json` returned `{"problems": []}`.
Editor diagnostics reported no errors in the four new gate/test files.

```sh
cd /home/ammslab3/softwarex_work/imq-ra-recovery
/home/ammslab3/softwarex_work/.venv/bin/python -m pytest -q tests/integration/test_clean_install_gate.py /home/ammslab3/softwarex_work/imq-umat-recovery/tests/test_clean_install_gate.py --import-mode=importlib
env PYTHONPATH=/home/ammslab3/softwarex_work/imq-ra-recovery:/home/ammslab3/softwarex_work/imq-umat-recovery/src:/home/ammslab3/otilib/build_py311 UMAT_OTI_REPO=/home/ammslab3/softwarex_work/imq-umat-recovery PYOTI_PATH=/home/ammslab3/otilib/build_py311 OTILIB_ROOT=/home/ammslab3/otilib/build_py311 RUN_OTILIB_TESTS=1 /home/ammslab3/softwarex_work/.venv/bin/python -m pytest -q -ra -m 'not abaqus and not arc and not network' --junitxml=.pytest_cache/recovery_install_offline.xml
cd /home/ammslab3/softwarex_work/imq-umat-recovery
env PYTHONPATH=/home/ammslab3/softwarex_work/imq-umat-recovery/src /home/ammslab3/softwarex_work/.venv/bin/python -m pytest -q -ra -m 'not abaqus and not arc and not network and not corpus_pass' --junitxml=.pytest_cache/recovery_install_offline.xml
```

## Remaining Failures And Prerequisites

1. Highest-priority real failure: UMAT
   `test_contract_fixtures.py::test_the_recorded_generation_is_this_worktrees_actual_transform`.
   Recorded transform fingerprint `94a92c01814f107a` differs from current
   `6aa20d22e37f14c9`. Updating a number alone would disguise stale evidence;
   shared generation/lock changes need genuine re-freezing or explicit retirement
   of affected frozen artifacts. Not changed in this packaging slice.
2. `test_the_two_denominators_stay_apart.py::test_every_terminal_state_in_the_report_is_marked_external_or_internal`:
   `arguments_diverged_before_the_routine` is EXTERNAL in the retained report
   but INTERNAL in the registry. Existing evidence needs regeneration/reconciliation.
3. External prerequisites: genuine matching ODB, licensed Abaqus 2021 odbAccess,
   compatible gfortran/Fortran runtime (9.4.0 here), and healthy Python. No private
   material source is required by the collaborator. No new analysis was launched.
4. Clean-clone release gate awaits an explicitly authorized commit and clean
   checkout of both final revisions. Generic UMAT/FCC, full-size models, OTILib
   binary redistribution, paper/corpus reproduction, other operating systems and
   the complete ledger remain outside this gate's verified scope.