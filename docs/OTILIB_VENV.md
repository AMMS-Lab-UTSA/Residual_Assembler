# Building OTILib for Python 3.11 without Conda

This page is the procedure for building a genuine OTILib (the hypercomplex
engine behind the OTI sensitivity paths) into a plain virtual environment,
without Conda. It is for users who need the OTILib-dependent workflows and
tests (the direct Python residual path, the finite-strain example, and the
tests run with `RUN_OTILIB_TESTS=1`) on a machine without Conda.
`scripts/setup_otilib.sh` is the Conda-based alternative; it stops when Conda
is unavailable. The procedure below was executed successfully during a
[clean-clone verification](evidence/recovery_clean_clone.md).

Prerequisites: Linux, Git, GCC/gfortran, GNU make, network access and a healthy
Python 3.11 with ctypes, ssl and venv. OTILib is an explicit external GPLv3
dependency, not the unrelated PyPI `pyoti` package. Do not reuse another
environment's extension modules or direction tables.

Use a fresh external work directory. `PY` must name the Python of the virtual
environment where Residual_Assembler is installed; `WORK` must not contain
either project checkout.

```sh
export WORK=/path/to/new/otilib_work      # a new directory outside both checkouts; the build stays here
export PY=/path/to/venv/bin/python        # the environment Residual_Assembler is installed in
mkdir "$WORK"
mkdir "$WORK/home"
export HOME="$WORK/home"
unset PYTHONPATH PYTHONHOME PYOTI_PATH OTILIB_ROOT
export PYTHONNOUSERSITE=1
export PIP_CONFIG_FILE=/dev/null
export PATH="$(dirname "$PY"):/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
"$PY" -I -m pip install cmake==3.31.10 Cython==3.3.0 numpy==2.4.6 scipy==1.17.1
git clone https://github.com/mauriaristi/otilib.git "$WORK/source"
git -C "$WORK/source" checkout --detach a4b7a05ca275e8d441b0b717b7b272e96728ffcc
git -C "$WORK/source" rev-parse HEAD
cmake -S "$WORK/source" -B "$WORK/build" -DBUILD_TESTING=OFF
cmake --build "$WORK/build" --target oticython --parallel 2
cmake --build "$WORK/build" --target gendata --parallel 2
export PYOTI_PATH="$WORK/build"
export OTILIB_ROOT="$WORK/build"
"$PY" -I -c 'from residual_core.algebra.otilib_adapter import otilib_status; status=otilib_status(); print(status); assert status["available"]; import pyoti.sparse as oti; value=(3+oti.e(1,order=2))**2; assert abs(value.get_deriv(1)-6)<1e-12; assert abs(value.get_deriv([1,1])-2)<1e-12; print(oti.__file__)'
```

Both selectors deliberately point to the **build directory**, not the source
root. The installed Residual_Assembler adapter adds that explicitly declared
external build to its import search path; no project-source `PYTHONPATH` is
required. Keep `PYOTI_PATH` and `OTILIB_ROOT` set to this build directory in
every shell that runs the OTILib workflows (the `HOME` override above is only
for the isolated build). CMake's
upstream `oticython` target uses the `python` executable on PATH and internally
compiles multiple extensions concurrently. The commands above run sequentially.
`BUILD_TESTING=OFF` disables OTILib's own upstream test suite; it does not skip
the project offline suites or the genuine sparse numerical probe.

The validated build used a fresh upstream clone and generated every extension
and direction table locally. GCC emitted warnings, but configure, build,
data generation and first/second derivative assertions all exited zero.
The exact toolchain, commands, dependencies and limitations are in the linked
evidence. This procedure does not establish support for other platforms or
other upstream revisions.