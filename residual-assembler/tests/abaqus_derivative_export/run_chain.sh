#!/usr/bin/env bash
# M3 full chain (real Abaqus): ODB -> derivative_fields.json -> resasm run -> du/da.
#
#   bash run_chain.sh [workdir]
#
# Reproduces the committed fixtures. Needs Abaqus (solver + python) and ifort on
# PATH. The Abaqus 2021 wrap-up may print a benign "buffer overflow detected"
# AFTER writing the ODB (a glibc-fortify abort in SMASimUtility); this script
# checks the .sta / .odb, not the exit code.
#
# IMPORTANT: the workdir must NOT contain parentheses or other shell-special
# characters -- Abaqus's UMAT compile step shells out and breaks on them (this
# repo's own path, ".../Residual_Assembler(2)/...", triggers it). The default
# workdir is therefore a paren-free temp dir, not a path under the repo.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
WORK="${1:-$(mktemp -d "${TMPDIR:-/tmp}/resasm_m3_XXXXXX")}"
UMAT="${2:-elastic_export_umat.for}"          # or elastic_oti_umat.for / cubic_oti_umat.for
INP="${3:-nonuniform_c3d8.inp}"               # or fcc_cubic_c3d8.inp
LAYOUT="${4:-derivative_layout.json}"         # or fcc_layout.json
ABQ="${ABAQUS_CMD:-abaqus}"
OTI="${OTI_DIR:-$HOME/MultiZ_f/oti}"

rm -rf "$WORK"; mkdir -p "$WORK"
cp "$HERE/$INP" "$WORK/model.inp"
cp "$HERE/$UMAT" "$HERE/$LAYOUT" "$WORK/"
# the exporter reads --layout by name; normalize to derivative_layout.json in WORK
[ "$LAYOUT" != "derivative_layout.json" ] && cp "$HERE/$LAYOUT" "$WORK/derivative_layout.json"
cp "$ROOT/residual_core/io/export_derivative_fields.py" \
   "$ROOT/residual_core/io/abaqus_odb_export.py" "$WORK/"
cd "$WORK" || exit 3

# The OTI UMAT (M4) links the MultiZ_f OTIM4N1 library via a local abaqus_v6.env.
case "$UMAT" in
  *oti*)
    if [ ! -f "$OTI/libotim4n1.a" ]; then
      echo "ERROR: OTI library not found at $OTI (set OTI_DIR)"; exit 2
    fi
    cp "$OTI/otim4n1.mod" "$OTI/libotim4n1.a" "$WORK/"
    cat > abaqus_v6.env <<PYENV
import os
oti = '$OTI'
compile_fortran += ['-I'+os.getcwd(), '-I'+oti]
link_sl  += ['-L'+oti, '-lotim4n1']
link_exe += ['-L'+oti, '-lotim4n1']
PYENV
    ;;
esac

echo "[1/3] Abaqus solve (C3D8 + $UMAT)"
"$ABQ" job=model user="$UMAT" interactive >/dev/null 2>&1
grep -q "COMPLETED SUCCESSFULLY" model.sta || { echo "ERROR: analysis did not complete"; exit 1; }

echo "[2/3] export ODB -> resasm_derivative_fields_v1"
"$ABQ" python export_derivative_fields.py -- \
    --odb model.odb --step Step-1 --frame -1 \
    --layout derivative_layout.json --output derivative_fields.json || exit 1

echo "[3/3] resasm run"
# parameters come from the layout metadata (single source of truth)
PARAMS=$(python3 -c "import json; print('[' + ', '.join(json.load(open('derivative_layout.json'))['sdv_layout']['parameters']) + ']')")
cat > sensitivity.yaml <<YAML
analysis:
  type: field_residual_sensitivity
  kinematics: small_strain
mesh:
  format: abaqus
  file: model.inp
results:
  format: resasm_derivative_fields_v1
  file: derivative_fields.json
parameters: $PARAMS
assumptions:
  parameter_independent_geometry: true
  parameter_independent_loads: true
  parameter_independent_boundaries: true
output:
  directory: results
YAML
python3 -c "import sys; sys.path.insert(0, '$ROOT'); from resasm_user.runner import run_from_config; r=run_from_config('sensitivity.yaml'); print('du/da norms:', r.summary['du_da_norms'])" || exit 1
echo "OK -> $WORK/results/"
