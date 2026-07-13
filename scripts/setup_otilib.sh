#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# setup_otilib.sh — install the genuine OTILib (GPLv3) as an EXTERNAL, optional
# dependency for the `--backend otilib` sensitivity path.
#
#   Repo:    https://github.com/mauriaristi/otilib.git   (branch master)
#   License: GPLv3  (external; NOT vendored into this framework's core)
#
# This is an explicit, opt-in helper. It is NEVER run by the test suite.
# Do NOT `pip install pyoti` — that PyPI name is an unrelated squat.
#
# Windows: OTILib builds only under WSL (per its README). Run this inside WSL.
#
# Usage:
#   scripts/setup_otilib.sh [TARGET_DIR]
#     TARGET_DIR  where to clone OTILib.
#                 default: <repo>/external/gpl/otilib   (clearly-marked GPL area)
#
# After it finishes, either `conda activate pyoti` in this shell, or export
#   OTILIB_ROOT=<TARGET_DIR>       (so the adapter can find pyoti)
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_URL="https://github.com/mauriaristi/otilib.git"
BRANCH="master"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FRAMEWORK_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
TARGET_DIR="${1:-${FRAMEWORK_ROOT}/external/gpl/otilib}"

echo "=============================================================="
echo " OTILib setup (GPLv3 external dependency)"
echo "   source : ${REPO_URL} (branch ${BRANCH})"
echo "   target : ${TARGET_DIR}"
echo "   NOTE   : OTILib is GPLv3. It is installed OUTSIDE this framework's"
echo "            core and only used when you pass --backend otilib."
echo "=============================================================="

# --- sanity checks --------------------------------------------------------
command -v git   >/dev/null 2>&1 || { echo "ERROR: git not found";   exit 1; }
command -v cmake >/dev/null 2>&1 || { echo "ERROR: cmake not found (install it or 'conda install cmake')"; exit 1; }
command -v conda >/dev/null 2>&1 || { echo "ERROR: conda not found (install Miniconda/Anaconda)"; exit 1; }

case "$(uname -s)" in
  Linux|Darwin) : ;;
  *) echo "WARNING: OTILib builds only on Unix/macOS/WSL. On Windows use WSL." ;;
esac

# --- 1. clone -------------------------------------------------------------
if [ -d "${TARGET_DIR}/.git" ]; then
  echo "[1/5] OTILib already cloned; pulling latest on ${BRANCH}."
  git -C "${TARGET_DIR}" fetch origin "${BRANCH}"
  git -C "${TARGET_DIR}" checkout "${BRANCH}"
  git -C "${TARGET_DIR}" pull --ff-only origin "${BRANCH}"
else
  echo "[1/5] Cloning OTILib..."
  mkdir -p "$(dirname "${TARGET_DIR}")"
  git clone --branch "${BRANCH}" "${REPO_URL}" "${TARGET_DIR}"
fi

cd "${TARGET_DIR}"

# --- 2. conda environment -------------------------------------------------
echo "[2/5] Creating/using conda environment 'pyoti'..."
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
if conda env list | grep -qE '^pyoti\s'; then
  echo "      conda env 'pyoti' already exists."
else
  conda env create -f environment.yml
fi
conda activate pyoti

# --- 3. configure + build -------------------------------------------------
echo "[3/5] Configuring + building with CMake..."
mkdir -p build
cd build
cmake ..
make

# --- 4. precomputed data --------------------------------------------------
echo "[4/5] Generating precomputed data (make gendata)..."
make gendata

# --- 5. put pyoti on the conda path --------------------------------------
echo "[5/5] Registering package on the conda path (conda develop .)..."
cd "${TARGET_DIR}"
conda develop . || conda develop "${TARGET_DIR}"

# --- report detected paths ------------------------------------------------
echo "=============================================================="
echo " OTILib build complete. Detected paths:"
PYPATH="$(python -c 'import pyoti, os; print(os.path.dirname(pyoti.__file__))' 2>/dev/null || true)"
if [ -n "${PYPATH}" ]; then
  echo "   python import path : ${PYPATH}"
  python -c "import pyoti.sparse as oti; x=3.5+oti.e(1,order=2); print('   self-test          : f.get_im(1) works ->', hasattr(x,'get_im'))"
else
  echo "   python import path : (pyoti not importable yet — run 'conda activate pyoti')"
fi
echo "   include dir        : ${TARGET_DIR}/include"
echo "   library/build dir  : ${TARGET_DIR}/build"
echo ""
echo " To use with this framework in a new shell, either:"
echo "     conda activate pyoti"
echo "   or:"
echo "     export OTILIB_ROOT=${TARGET_DIR}"
echo "     export OTILIB_INCLUDE_DIR=${TARGET_DIR}/include"
echo "     export OTILIB_LIBRARY_DIR=${TARGET_DIR}/build"
echo ""
echo " Verify:"
echo "   python -c \"from residual_core.algebra.otilib_adapter import otilib_status; print(otilib_status())\""
echo "=============================================================="
