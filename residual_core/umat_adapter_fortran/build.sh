#!/usr/bin/env bash
#=======================================================================
#  build.sh  -- gfortran build recipe for the standalone UMAT replay
#               adapter.  Honest about what compiles and what needs
#               ifort+Abaqus.  Runs on Git-Bash/MSYS2 or Linux.
#
#  Produces (in ./build):
#    umat_driver_mock   -- driver + mock UMAT.  ALWAYS builds with
#                          gfortran.  Use it to test the plumbing and
#                          umat_replay.py --dry-run.
#    umat_driver        -- driver + REAL Grilli UMAT + stubs.  Only
#                          links if umat.for compiles (needs ifort+Abaqus
#                          or a patched source tree -- see README.md).
#
#  Usage:
#    ./build.sh [path-to-ngrilli_Oxford_Crystal_Plasticity]
#  Default source path is ../../sources/permissive/ngrilli_Oxford_Crystal_Plasticity
#=======================================================================
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="${1:-$HERE/../../sources/permissive/ngrilli_Oxford_Crystal_Plasticity}"
STUBS="$HERE/aba_stubs"
BUILD="$HERE/build"
FC="${FC:-gfortran}"
EXE=""
case "$(uname -s)" in MINGW*|MSYS*|CYGWIN*) EXE=".exe";; esac

mkdir -p "$BUILD"
cd "$BUILD"

echo "== toolchain =="
$FC --version | head -1
echo "SRC   = $SRC"
echo "STUBS = $STUBS"
echo

if [ ! -f "$SRC/umat.for" ]; then
  echo "WARNING: $SRC/umat.for not found; only the mock path will build."
fi

# common flags
FREE="-O2 -fbacktrace"
FIXED="-O2 -ffixed-line-length-none"
# flags to give the unmodified ifort-oriented UMAT its best chance on gfortran
UMATFLAGS="-c -cpp -fcray-pointer -ffixed-line-length-none -fno-range-check \
 -std=legacy -fallow-argument-mismatch -fallow-invalid-boz -w \
 -I$STUBS -I$SRC"

set -x
# ---- 1. driver + mock + stubs (our code: must all compile) ----
$FC $FREE  -c "$HERE/umat_driver.f90"      -o umat_driver.o
$FC $FREE  -c "$STUBS/umat_mock.f90"       -o umat_mock.o
$FC $FIXED -c "$STUBS/aba_stubs.f"         -o aba_stubs.o
$FC $FIXED -c "$STUBS/lapack_stub.f"       -o lapack_stub.o
set +x

echo
echo "== link umat_driver_mock (plumbing test binary) =="
if $FC $FREE umat_driver.o umat_mock.o -o "umat_driver_mock$EXE"; then
  echo "OK  -> build/umat_driver_mock$EXE"
else
  echo "FAILED to link umat_driver_mock (unexpected -- these are our files)."
fi

echo
echo "== attempt to compile the REAL Grilli UMAT (umat.for) =="
echo "   log -> build/umat_compile.log"
if [ -f "$SRC/umat.for" ]; then
  if $FC $UMATFLAGS "$SRC/umat.for" -o umat.o 2> umat_compile.log; then
    echo "umat.for COMPILED. Linking real umat_driver ..."
    if $FC $FREE umat_driver.o umat.o aba_stubs.o lapack_stub.o \
         -o "umat_driver$EXE"; then
      echo "OK  -> build/umat_driver$EXE  (REAL crystal-plasticity UMAT)"
    else
      echo "umat.o built but link failed -- see missing symbols above."
    fi
  else
    echo "umat.for did NOT compile with gfortran (expected)."
    echo "First errors:"
    grep -E "Error:|Fatal" umat_compile.log | head -12 | sed 's/^/   /'
    echo "   ... full log in build/umat_compile.log"
    echo
    echo "This is the documented gfortran limitation (see README.md):"
    echo "  * Cray POINTEE + TARGET on the twin arrays (ifort extension)"
    echo "  * REAL(4)/REAL(8) mismatch on function 'trace'"
    echo "  Build the real UMAT with Intel ifort + Abaqus (MKL), then link:"
    echo "    ifort ... -c umat.for   (via 'abaqus make' / the Abaqus env)"
    echo "    \$FC umat_driver.o umat.o aba_stubs.o -o umat_driver   (or MKL)"
  fi
else
  echo "SKIP: umat.for not present at \$SRC."
fi

echo
echo "== done. binaries in: $BUILD =="
ls -1 "$BUILD"/umat_driver* 2>/dev/null || true
