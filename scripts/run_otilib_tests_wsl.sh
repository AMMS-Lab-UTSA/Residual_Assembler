#!/usr/bin/env bash
#=============================================================================
# run_otilib_tests_wsl.sh — deterministically activate the genuine OTILib
# (GPLv3, external) and run the three OTILib test files.
#
# ONE copy-paste command. It either RUNS the real OTILib tests, or exits
# non-zero and tells you EXACTLY what is missing and how to fix it. It never
# lets the OTILib tests silently skip: it exports RUN_OTILIB_TESTS=1, which
# turns the skip guard into a hard assert, so a green run here is real numeric
# validation and nothing else can be mistaken for one.
#
#   From inside WSL/Linux:   bash scripts/run_otilib_tests_wsl.sh
#   From Windows:            powershell -File scripts\run_otilib_tests_wsl.ps1
#
# Overridable (all auto-detected if unset):
#   OTILIB_ROOT         path of the OTILib source build   (e.g. /root/otilib)
#   OTILIB_CONDA_ENV    conda env holding the built pyoti (default: pyoti)
#   CONDA_SH            path to conda's profile.d/conda.sh
#
# OTILib is NOT vendored here and is NOT a dependency of this framework. It is
# GPLv3 and must be installed separately: scripts/setup_otilib.sh
# (https://github.com/mauriaristi/otilib.git). Do NOT `pip install pyoti` —
# that PyPI name is an unrelated package.
#=============================================================================
set -u

TESTS=(
  "tests/framework/test_otilib_adapter.py"
  "tests/framework/test_otilib_spring_sensitivity.py"
  "tests/framework/test_otilib_fe_sensitivity.py"
)

die() {
  echo ""
  echo "CANNOT RUN THE OTILIB TESTS"
  echo "--------------------------"
  echo "$1"
  echo ""
  echo "OTILib is an optional, external, GPLv3 library. Install it with:"
  echo "    bash scripts/setup_otilib.sh"
  echo "Then re-run this script. (Do NOT 'pip install pyoti' — unrelated package.)"
  exit 1
}

# --- repo root: derived from this script's own location (works at /mnt/c/... too)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$HERE")"
cd "$REPO" || die "cannot cd into the repo root: $REPO"

echo "== OTILib test run (deterministic activation) =="
echo "repo            : $REPO"

# --- 1. conda ---------------------------------------------------------------
if [ -z "${CONDA_SH:-}" ]; then
  for c in "${CONDA_PREFIX:-}/etc/profile.d/conda.sh" \
           "$HOME/miniconda3/etc/profile.d/conda.sh" \
           "$HOME/anaconda3/etc/profile.d/conda.sh" \
           "/root/miniconda3/etc/profile.d/conda.sh" \
           "/opt/conda/etc/profile.d/conda.sh"; do
    [ -f "$c" ] && { CONDA_SH="$c"; break; }
  done
fi
if [ -z "${CONDA_SH:-}" ] && command -v conda >/dev/null 2>&1; then
  CONDA_SH="$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh"
fi
[ -n "${CONDA_SH:-}" ] && [ -f "$CONDA_SH" ] \
  || die "conda not found. Looked for profile.d/conda.sh under \$CONDA_PREFIX,
\$HOME/miniconda3, \$HOME/anaconda3, /root/miniconda3, /opt/conda, and on PATH.
Set CONDA_SH=/path/to/conda.sh if it lives somewhere else."
# shellcheck disable=SC1090
source "$CONDA_SH"
echo "conda           : $CONDA_SH"

# --- 2. env holding the built pyoti ----------------------------------------
ENV_NAME="${OTILIB_CONDA_ENV:-pyoti}"
conda activate "$ENV_NAME" 2>/dev/null \
  || die "conda env '$ENV_NAME' not found (conda env list).
Create it with scripts/setup_otilib.sh, or set OTILIB_CONDA_ENV=<name>."
echo "conda env       : $ENV_NAME  ($(python --version 2>&1))"

# --- 3. OTILIB_ROOT ---------------------------------------------------------
if [ -z "${OTILIB_ROOT:-}" ]; then
  for r in "/root/otilib" "$HOME/otilib" "$REPO/external/gpl/otilib" "/opt/otilib"; do
    [ -d "$r" ] && { OTILIB_ROOT="$r"; break; }
  done
fi
[ -n "${OTILIB_ROOT:-}" ] && [ -d "$OTILIB_ROOT" ] \
  || die "OTILib source build not found. Looked in /root/otilib, \$HOME/otilib,
$REPO/external/gpl/otilib, /opt/otilib. Set OTILIB_ROOT=/path/to/otilib."
export OTILIB_ROOT
echo "OTILIB_ROOT     : $OTILIB_ROOT"

# The adapter probes \$OTILIB_ROOT, \$OTILIB_ROOT/src/python and \$OTILIB_ROOT/build,
# so OTILIB_ROOT alone is enough. Add build/ to PYTHONPATH anyway for anything
# that imports pyoti without going through the adapter.
export PYTHONPATH="${OTILIB_ROOT}/build${PYTHONPATH:+:$PYTHONPATH}"

# --- 4. the library must actually import (and be the GENUINE OTI) -----------
python - <<'PY' || die "the genuine OTILib (pyoti) is present on disk but does not import
in this env. Rebuild it: bash scripts/setup_otilib.sh"
import sys
sys.path.insert(0, ".")
from residual_core.algebra.otilib_adapter import otilib_available, otilib_status
st = otilib_status()
if not otilib_available():
    print("otilib_available() = False:", st.get("error", "")[:200])
    raise SystemExit(1)
print("pyoti module    : %s (GENUINE OTI API verified)" % st.get("api_module"))
PY

# --- 5. pytest --------------------------------------------------------------
python -c "import pytest" 2>/dev/null \
  || die "pytest is not installed in conda env '$ENV_NAME'. Install it with:
    conda run -n $ENV_NAME python -m pip install pytest"

# --- 6. run the real thing --------------------------------------------------
# RUN_OTILIB_TESTS=1 => the OTILib tests must RUN; if OTILib were missing they
# would FAIL rather than skip. A green result here is genuine numeric validation.
export RUN_OTILIB_TESTS=1
echo "RUN_OTILIB_TESTS: 1  (skip guard disabled -- these must really run)"
echo ""
python -m pytest "${TESTS[@]}" -q -rs
rc=$?
echo ""
if [ $rc -eq 0 ]; then
  echo "RESULT: OTILib tests PASSED against the genuine OTILib at $OTILIB_ROOT"
else
  echo "RESULT: OTILib tests FAILED (exit $rc) -- this is a real failure, not a skip."
fi
exit $rc
