# Recovery W1 evidence (2026-09-18)

## Scope and status

Checkout: `/home/ammslab3/softwarex_work/imq-ra-recovery`.
Branch: `integration/imqcam-recovery-2026-09-18`.
Read the shared `imq_BRIEF.md` and `docs/BRANCH_IMPLEMENTATION_AUDIT.md`.
No agents, commits, pushes, branch changes, Abaqus jobs, network jobs, or
whole-corpus runs. Original Claude worktrees were not changed. Other workers'
changes in `resasm_user/__init__.py`, `residual_core/ui/cli.py`, and
`tests/framework/test_user_layer.py` were preserved; none was edited by W1.
No replay, IO, UI, report implementation, or packaging files were edited.

**Clean-install verification: NOT RUN / NOT ESTABLISHED.** This is validation
against an existing external OTILib build, not evidence of a reproducible clean
installation. No completion-ledger rows are claimed complete.

Bounded W1 result: both documented baseline failures resolved; all 15 OTILib
tests passed without skips; full RA offline suite completed with 296 passed,
19 skipped, and zero failures/errors. There are no remaining observed failures
in the selected offline suite; skips and out-of-scope capabilities remain below.

## Interpreter and genuine OTILib

All Python commands use `/home/ammslab3/softwarex_work/.venv/bin/python`
(Python 3.11.7), from the recovery cwd, with explicit import paths:

```sh
cd /home/ammslab3/softwarex_work/imq-ra-recovery
export PYTHONPATH=/home/ammslab3/softwarex_work/imq-ra-recovery:/home/ammslab3/otilib/build_py311
export OTILIB_ROOT=/home/ammslab3/otilib/build_py311
export RUN_OTILIB_TESTS=1
```

The exports above abbreviate the identical inline environment assignments used
in the validation commands. The initial import/baseline probe preceded setting
`RUN_OTILIB_TESTS=1`; all post-edit validation required it.

Before using the build, ran these read-only checks:

```sh
git -C /home/ammslab3/otilib status --short --branch
git -C /home/ammslab3/otilib branch --show-current
ls -ld /home/ammslab3/otilib/build_py311
rg --files /home/ammslab3/otilib/build_py311 -g '*.so' -g 'pyoti*' -g '*config*' -g '*log*'
```

Observed a clean OTILib checkout on `py311-numpy2`, with CPython-3.11 extension
modules for `core`, `sparse`, `dense`, `real`, and `fem`. Reused that build
without modifying or rebuilding the external clone. CMake's existing cache
records GCC 9.4, `GENERATE_CYTHON_OTI=ON`, `GENERATE_FORTRAN=OFF`, and
`USE_OPENMP=OFF`. The handoff attributes the build to Cython 0.29.37; no claim
is made that W1 reran those build commands.

Actual import/derivative probe:

```sh
/home/ammslab3/softwarex_work/.venv/bin/python -c 'import sys,residual_core,pyoti,pyoti.sparse; from residual_core.algebra.otilib_adapter import otilib_status,OtiContext; print(sys.executable); print(residual_core.__file__); print(pyoti.__file__); print(pyoti.sparse.__file__); print(otilib_status()); context=OtiContext(2,3); value=context.seed(2.,1)**3; print("cubic third derivative:",context.deriv(value,(3,0))); assert context.deriv(value,(3,0)) == 6.'
```

Observed:

```text
/home/ammslab3/softwarex_work/.venv/bin/python
/home/ammslab3/softwarex_work/imq-ra-recovery/residual_core/__init__.py
/home/ammslab3/otilib/build_py311/pyoti/__init__.py
/home/ammslab3/otilib/build_py311/pyoti/sparse.cpython-311-x86_64-linux-gnu.so
available=True; api_module=pyoti.sparse; error=''
cubic third derivative: 6.0
```

This is the genuine compiled OTI API, not the unrelated PyPI package. No
`pip install pyoti`, fallback backend, or dependency installation was used.

## Baseline fixes and history

Inspected `git log` and `git show ab0ee10 --
residual_core/core/where_it_went_wrong.py` before changing expectations. That WIP
commit added two real diagnostic layers and stopped treating a secant across
an evolving-state transition as proof of a bad tangent. The two failing tests
still asserted the superseded behavior. Production diagnosis was not weakened.

- Layer regression now pins all eleven names, their order, and owner coverage,
  including `parameter_sensitivity` and `element_jacobian`.
- J2 regression now requires `constitutive_derivative=not_established`, positive
  state movement, a measured secant disagreement greater than 0.1, an explicit
  request for `material_update`, and `complete=False`. UMAT and transformation
  evidence must still hold, and downstream layers must still be checked.
- The first-yield regression retains the measured increment check and now also
  requires the constitutive result to remain unestablished without an update.

Initial reproduction command:

```sh
/home/ammslab3/softwarex_work/.venv/bin/python -m pytest -q -ra tests/framework/test_otilib* tests/framework/test_a_verified_umat_drives_the_assembler.py::test_a_verified_fixture_passes_every_stage tests/framework/test_which_layer_the_failure_came_from.py::test_the_nine_layers_are_the_nine_the_brief_names
```

Result: **2 failed, 6 passed, 0 skipped**. Failures were exactly the documented
J2 verdict and nine-versus-eleven layer expectations. The two tests have been
renamed to describe the corrected contracts.

## Scalar-generic FE propagation

Inspected the local C3D8 history and
`origin/cross-platform-hardening:residual-assembler/residual_core/formulations/solid3d_kernel.py`
before editing. No alternative element family or physics implementation was
imported. The first new analytic force tests reproduced two failures at the
`dtype=float` stress conversion, then passed after the minimal dtype change.

W6 coordination note: only scalar-preserving array construction and accumulator
promotion changed in the existing C3D8 kernels. Geometry, quadrature, tensor
conventions, formulas, and public signatures are unchanged.

- `residual_core/formulations/c3d8_kernel.py`: retain stress and material-tangent
  scalar types through small/finite-strain force, Voigt mapping, material and
  geometric stiffness, and global force scatter.
- `residual_core/formulations/c3d8_sensitivity.py`: retain derivative scalar
  types through element `dR/dp` and global scatter.
- `tests/framework/test_otilib_fe_sensitivity.py`: nine added cases cover small
  and finite strain, mixed real/OTI contributions, nonsequential node IDs,
  material stiffness, and geometric stiffness with real D and OTI stress.

Independent reference: stress amplitude `p**2*q`, seeded at `p=2`, `q=3`, with
two OTI bases and order three. Its real value and derivatives for directions
`(0,0), (1,0), (0,1), (2,0), (1,1), (2,1)` are respectively
`12, 12, 4, 6, 4, 2`. Nodal forces are checked against analytic cube-face
tractions (including stretched face areas), not another call to the volume
kernel. Stiffness acts on an affine field with an independently specified
constant stress. Both relative and absolute tolerances are `1e-13`.

Runnable API example/verification with inputs, expected derivatives, and
direction-specific mismatch diagnostics is the new C3D8 test group:

```sh
/home/ammslab3/softwarex_work/.venv/bin/python -m pytest -q -ra tests/framework/test_otilib_fe_sensitivity.py -k c3d8
```

The existing spring and nonlinear two-element bar tests additionally compare
OTILib sensitivities against independently re-solved finite differences. No
new GUI/CLI command was added: UI and registration are owned by other workers.

## Validation commands and counts

All counts below separate skips from passes. Terminal output was occasionally
replaced by unrelated concurrent-worker output; those results were not counted
as W1 evidence. Later suites use W1-specific JUnit files under `.pytest_cache`.

| Check | Passed | Failed/errors | Skipped | Notes |
| --- | ---: | ---: | ---: | --- |
| Initial baseline plus existing OTILib tests | 6 | 2 | 0 | Expected reproduction |
| Baseline correction plus existing OTILib tests | 40 | 0 | 0 | Before FE changes |
| Initial new force regression probe | 0 | 2 | 0 | 2 deselected; float coercion |
| Force regression after fix | 2 | 0 | 0 | 2 deselected |
| FE propagation plus baseline and real C3D8 | 52 | 0 | 0 | Before geometric-stiffness case |
| Final focused tests | 57 | 0 | 0 | All OTILib files included |
| Framework offline | 214 | 0 | 8 | 222 executed; JUnit verified |
| Full offline | 296 | 0 | 19 | 315 executed; final JUnit verified |
| OTILib subset of full offline | 15 | 0 | 0 | Every `test_otilib_*` case executed |

Final focused command:

```sh
/home/ammslab3/softwarex_work/.venv/bin/python -m pytest -q -ra tests/framework/test_otilib_adapter.py tests/framework/test_otilib_fe_sensitivity.py tests/framework/test_otilib_spring_sensitivity.py tests/framework/test_a_verified_umat_drives_the_assembler.py tests/framework/test_which_layer_the_failure_came_from.py tests/framework/test_c3d8_sensitivity.py
```

Framework command (actual Python wrapper asserts cwd and prints import origin):

```sh
/home/ammslab3/softwarex_work/.venv/bin/python -c 'import os,sys,pytest,residual_core; assert os.getcwd()=="/home/ammslab3/softwarex_work/imq-ra-recovery",os.getcwd(); print("W1 FRAMEWORK",os.getcwd(),sys.executable,residual_core.__file__,flush=True); sys.exit(pytest.main(["-q","-ra","tests/framework","-m","not abaqus and not arc and not network","--junitxml=.pytest_cache/recovery_W1_framework.xml"]))'
```

`recovery_W1_framework.xml`: 222 tests, 8 skips, 0 failures, 0 errors,
18.861 seconds. Editor diagnostics: no errors in the five changed Python files.

Full-suite attempts wrote `recovery_W1_offline.xml` (187 executed, 18 skipped)
and `recovery_W1_offline_complete.xml` (123 executed, 11 skipped), but both
were interrupted with status 130. Neither is completion evidence.

The terminal was shared with other workers and returned unrelated output or
SIGINT while commands ran. The full run was therefore isolated from the
terminal's process group. Final actual command (environment shown explicitly):

```sh
cd /home/ammslab3/softwarex_work/imq-ra-recovery
setsid --wait env PYTHONPATH=/home/ammslab3/softwarex_work/imq-ra-recovery:/home/ammslab3/otilib/build_py311 OTILIB_ROOT=/home/ammslab3/otilib/build_py311 RUN_OTILIB_TESTS=1 /home/ammslab3/softwarex_work/.venv/bin/python -m pytest -q -ra -m 'not abaqus and not arc and not network' --junitxml=/home/ammslab3/softwarex_work/imq-ra-recovery/.pytest_cache/recovery_W1_offline_final.xml
```

Although the shared terminal returned status 130, the isolated pytest process
finished and wrote a complete JUnit report: 315 test cases, 0 failures, 0 errors,
19 skips, 19.208 seconds, timestamp `2026-09-18T09:19:35.605027-05:00`.
An earlier isolated run also completed (same counts, 18.714 seconds) in
`recovery_W1_offline_isolated.xml`, despite unrelated terminal output. A process
inspection confirmed no W1 pytest process remained. Completion is based on the
finished XML reports, not the misleading terminal status or a partial report.

Parsed the final report using `xml.etree.ElementTree`, checked every testcase
whose classname contains `.test_otilib_`, asserted exactly 15 cases and zero
skip/failure/error children, and counted all skip reasons. Machine-readable
summary of that measured result:

```json
{
  "executed": 315,
  "passed": 296,
  "failed": 0,
  "errors": 0,
  "skipped": 19,
  "otilib": {"passed": 15, "failed": 0, "skipped": 0},
  "skips": {"missing_frozen_store": 11, "uninitialized_external_source": 8},
  "clean_install_verified": false
}
```

The 11 frozen-store skips require
`corpus_run/pass11/results/store_verification.jsonl`; they explicitly mean the
two readers were not compared on those real records. The other 8 require
`sources/permissive/ngrilli_Oxford_Crystal_Plasticity/ExampleInputFiles/HCPnoTwin/Compression111.inp`.
Neither dependency was fetched or reconstructed, and no skip was counted as
a pass. Final read-only OTILib status still reports `py311-numpy2`, clean.

## Remaining boundaries

Geometry and displacement inputs still use real NumPy linear algebra. This is
generic stress/tangent propagation on fixed real geometry, not generic shape
sensitivities or a full arbitrary-order nonlinear FE solve. The explicit
`solve_du_dp` remains a real linear solve after coefficient extraction.
Other element families, UMAT provider/replay integration, and GUI wiring are
not claimed by this slice. J2 constitutive derivatives without a re-evaluated
update remain NOT ESTABLISHED by design. External-source skips are not passes.

## Files changed by W1

1. `residual_core/formulations/c3d8_kernel.py`
2. `residual_core/formulations/c3d8_sensitivity.py`
3. `tests/framework/test_a_verified_umat_drives_the_assembler.py`
4. `tests/framework/test_which_layer_the_failure_came_from.py`
5. `tests/framework/test_otilib_fe_sensitivity.py`
6. `docs/evidence/recovery_W1.md`