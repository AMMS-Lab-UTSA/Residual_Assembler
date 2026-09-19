#!/usr/bin/env bash
# Reproduce Residual_Assembler and UMAT_source_transformation from fresh clones
# of their published main branches, following only the documented
# instructions, and record every step.
#
# Usage (docs/INSTALL.md, section 7):
#
#   export PYOTI_PATH=/path/to/otilib/build OTILIB_ROOT=/path/to/otilib/build
#   bash scripts/reproduce_from_clean_clones.sh --python python3.11 \
#       --odb /path/to/Analysis.odb --cantilever /path/to/cantilever/work NEW_DIR
#
# Inputs, each made as the documents say:
#   PYOTI_PATH, OTILIB_ROOT  the OTILib build directory (docs/INSTALL.md section 4)
#   --python                 a Python 3.10 or newer with venv, ctypes and ssl (section 1)
#   --odb                    the Abaqus result of examples/presentation_request/Analysis.inp (section 7)
#   --cantilever             a folder holding j2/cantilever_j2_nominal.inp and .odb (section 7, item 8)
#   NEW_DIR                  a folder that does not exist yet
#
# The user steps, in order, each with its section of docs/INSTALL.md:
#   clone_ra, clone_umat  git clone of both repositories side by side (section 2)
#   gate                  scripts/clean_install_gate.py --branch main --cantilever (section 7)
#   suite_ra, suite_umat  both offline suites in the environment the gate made (sections 6 and 7)
#   examples              scripts/audit_recovery_usage.py --phase examples (section 6, item 7)
# Each step's exit code, time and log go to NEW_DIR/summary.tsv; the suites'
# JUnit XML to ra_suite.xml and umat_suite.xml; the gate's report to
# gate/report.json; the examples' records to usage/.
#
# Harness-only bookkeeping, which a user does not need, is marked "harness"
# below: the run starts again without the calling shell's Python settings, and
# `record` saves the output of a read-only query (the commits of the clones and
# of the OTILib build, and whether the run left the clones unmodified) and
# lists it in NEW_DIR/harness.tsv. docs/evidence/clean_clone_commands.json lists
# every step, variable and harness line of this script with the instruction it
# rests on (or, for a harness line, why a user does not need it);
# `python tools/audit_integrity.py --check ci3` checks that list against this file.
set -u

# harness: start again without the calling shell's Python settings (PYTHON*,
# PIP_*, CONDA*, VIRTUAL_ENV) and without the variables this script sets itself,
# so that another Python installation of the operator's cannot reach the run;
# everything else (PATH, compilers, licences) is kept
if [ -z "${REPRODUCE_CLEAN_SHELL:-}" ]; then
    exec env $(compgen -e | grep -E '^(PYTHON|PIP_|CONDA)|^(VIRTUAL_ENV|UMAT_OTI_REPO|RUN_OTILIB_TESTS)$' | sed 's/^/-u /') \
        REPRODUCE_CLEAN_SHELL=1 bash "$0" "$@"
fi

usage() { sed -n '6,10p' "$0" >&2; exit 2; }
PY=; ODB=; CANT=
while [ $# -gt 1 ]; do
    case $1 in
        --python) PY=$2; shift 2 ;;
        --odb) ODB=$(realpath "$2"); shift 2 ;;
        --cantilever) CANT=$(realpath "$2"); shift 2 ;;
        *) usage ;;
    esac
done
[ $# -eq 1 ] && [ -n "$PY" ] && [ -n "$ODB" ] && [ -n "$CANT" ] || usage
C=$(realpath -m "$1")
[ -e "$C" ] && { echo "refusing: $C exists; give a new folder" >&2; exit 2; }
[ -n "${PYOTI_PATH:-}" ] && [ -n "${OTILIB_ROOT:-}" ] || {
    echo "set PYOTI_PATH and OTILIB_ROOT to the OTILib build directory (docs/INSTALL.md section 4)" >&2; exit 2; }
mkdir -p "$C/logs"
RA=$C/Residual_Assembler
UMAT=$C/UMAT_source_transformation
printf "step\texit\tseconds\tlog\n" > "$C/summary.tsv"
printf "record\texit\tfile\n" > "$C/harness.tsv"

step() {  # step NAME CWD COMMAND...: a documented user step
    local name=$1 cwd=$2; shift 2
    local started=$SECONDS
    (cd "$cwd" && "$@") > "$C/logs/$name.log" 2>&1
    local status=$?
    printf "%s\t%s\t%s\t%s\n" "$name" "$status" "$((SECONDS - started))" "logs/$name.log" >> "$C/summary.tsv"
    echo "$name: exit $status"
    return $status
}

record() {  # record FILE COMMAND...: harness only, a read-only query saved to NEW_DIR/FILE
    local file=$1; shift
    "$@" > "$C/$file" 2>> "$C/logs/harness.log"
    local status=$?
    printf "%s\t%s\t%s\n" "$*" "$status" "$file" >> "$C/harness.tsv"
    return $status
}

finish() {
    record ra_status_after.txt git -C "$RA" status --porcelain --untracked-files=all
    record umat_status_after.txt git -C "$UMAT" status --porcelain --untracked-files=all
    if [ -s "$C/ra_status_after.txt" ] || [ -s "$C/umat_status_after.txt" ]; then
        echo "the run changed a clone: see ra_status_after.txt and umat_status_after.txt" >&2
    fi
    cat "$C/summary.tsv"
    awk -F'\t' 'NR > 1 && $2 != 0 { failed = 1 } END { exit failed }' "$C/summary.tsv"
    exit $?
}

# section 2: both repositories side by side (the default branch is main)
step clone_ra "$C" git clone https://github.com/AMMS-Lab-UTSA/Residual_Assembler.git || finish
step clone_umat "$C" git clone https://github.com/AMMS-Lab-UTSA/UMAT_source_transformation.git || finish
record ra_head.txt git -C "$RA" rev-parse HEAD
record umat_head.txt git -C "$UMAT" rev-parse HEAD

# harness: the OTILib build must be at the commit section 4 pins
record otilib_head.txt git -C "$PYOTI_PATH" rev-parse HEAD
record otilib_pin.txt grep -o -m 1 'checkout --detach [0-9a-f]\{40\}' "$RA/docs/INSTALL.md"
if [ "checkout --detach $(cat "$C/otilib_head.txt")" != "$(cat "$C/otilib_pin.txt")" ]; then
    echo "refusing: the OTILib build at $PYOTI_PATH is not at the commit docs/INSTALL.md section 4 pins" >&2
    exit 2
fi

# section 7: the clean-install gate, from the parent folder of the clones
step gate "$C" "$PY" Residual_Assembler/scripts/clean_install_gate.py \
    --umat-repo UMAT_source_transformation --python "$PY" \
    --odb "$ODB" --work "$C/gate" --branch main --cantilever "$CANT"
[ -f "$C/gate/env/bin/activate" ] || finish

# sections 6 and 7: the environment the gate made, with OTILib (section 4); the
# packages are installed from wheels, so the tests read the UMAT models from
# the checkout named by UMAT_OTI_REPO (section 6, item 6)
. "$C/gate/env/bin/activate"
export PYOTI_PATH OTILIB_ROOT
export RUN_OTILIB_TESTS=1
export UMAT_OTI_REPO="$UMAT"
step suite_ra "$RA" python -m pytest -q -m "not abaqus and not arc and not network" --junitxml="$C/ra_suite.xml"
step suite_umat "$UMAT" python -m pytest -q --junitxml="$C/umat_suite.xml"

# section 6, item 7: every example that needs no Abaqus, in one command
step examples "$RA" python scripts/audit_recovery_usage.py --umat ../UMAT_source_transformation \
    --phase examples --work "$C/examples" --evidence-dir "$C/usage"

finish
