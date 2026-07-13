#!/usr/bin/env bash
#
# demo_framework.sh — end-to-end demo of the model-agnostic residual assembler.
#
# Runs entirely offline (no Abaqus, no external code). Uses the module form of the
# CLI so it works whether or not the `resasm` console script is installed.
#
#   bash scripts/demo_framework.sh
#
set -euo pipefail

# repo root = parent of this script's directory
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# Interpreter is configurable so it works on any host: set PYTHON if `python`
# is not on PATH (e.g. PYTHON=python3, PYTHON=py, or a full path).
PYTHON="${PYTHON:-python}"
RESASM="$PYTHON -m residual_core.ui.cli"
EX="residual_core/examples"

hr() { printf '\n=== %s ===\n' "$1"; }

hr "1. Registered backends (crystal plasticity is ONE of them)"
$RESASM backends

hr "2. Regenerate the minimal example assets"
$PYTHON -m residual_core.examples.generate_minimal

hr "3. Inspect a truss model (auto-detection)"
$RESASM inspect "$EX/minimal_truss/model.json"

hr "4. Assemble the truss residual (Mode 4, formulation)"
$RESASM assemble "$EX/minimal_truss/model.json" --mode formulation

hr "5. Inspect a mixed truss+beam model with per-element detail"
$RESASM inspect "$EX/minimal_mixed/model.json" --detail

hr "6. Requirements for a solid stress-driven model with NO field yet"
$RESASM requirements "$EX/minimal_c3d8_stress_driven/model.json" --mode stress-driven || true

hr "7. Assemble the solid residual once the field is attached (Mode 1)"
$RESASM assemble "$EX/minimal_c3d8_stress_driven/model.json" \
        --mode stress-driven --fields "$EX/minimal_c3d8_stress_driven/fields.json"

hr "8. External residual comparison (offline; needs only the field export)"
$PYTHON scripts/compare_residuals.py \
        --model "$EX/minimal_c3d8_stress_driven/model.json" \
        --fields "$EX/minimal_c3d8_stress_driven/fields.json"

hr "9. Abaqus scripts skip cleanly when Abaqus is absent"
$PYTHON scripts/extract_odb_fields.py --odb missing.odb || true
$PYTHON scripts/compare_uel_rhs.py || true

hr "10. Offline framework test suite"
for t in test_assembler test_truss2_backend test_beam2_backend \
         test_mixed_model_dispatch test_dof_manager_mixed \
         test_requirements_negative test_public_api test_neutral_io; do
  $PYTHON "tests/framework/$t.py" >/dev/null && echo "PASS $t"
done

printf '\nDemo complete. See REVIEW_GUIDE.md for the guided tour.\n'
