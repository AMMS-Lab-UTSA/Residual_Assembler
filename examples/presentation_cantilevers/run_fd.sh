#!/usr/bin/env bash
# Perturbed-parameter Abaqus reruns for the central-difference reference.
#   run_fd.sh j2|fcc "<rel steps>"
# One Abaqus job at a time. Success is judged from the .sta line
# "THE ANALYSIS HAS COMPLETED SUCCESSFULLY": Abaqus 2021.HF5 aborts with
# signal 6 during teardown after writing a complete ODB, so its exit code is
# not the verdict. Each run is exported to <job>_fields.npz and logged in
# fd_manifest.tsv (job, parameter, sign, rel step, value, status, seconds).
#
# The UMAT is the ORIGINAL source from UMAT_source_transformation:
#   j2  -> parameter_sensitivity/models/m3_j2/umat.for
#   fcc -> parameter_sensitivity/models/m6_fcc/umat.for
# set UMAT=/path/to/umat.for, or copy it to <model>/<model>_umat.for.
# Runs happen in ./<model>/ below the current directory (WORK overrides);
# keep that outside any repository: each rerun ODB is exported and deleted.
set -u
model=$1; steps=${2:-"0.02 0.01 0.005"}
here=$(cd "$(dirname "$0")" && pwd)
work=${WORK:-$PWD}
mkdir -p "$work/$model" && cd "$work/$model"
umat=${UMAT:-${model}_umat.for}
[ -f "$umat" ] || { echo "UMAT source not found: $umat (set UMAT=...)" >&2; exit 2; }
read -r -a names < <(python3 -c "import sys; sys.path.insert(0,'$here'); from gen_cantilever import MODELS; print(' '.join(MODELS['$model']['prop_names']))")
read -r -a nominal < <(python3 -c "import sys; sys.path.insert(0,'$here'); from gen_cantilever import MODELS; print(' '.join(repr(p) for p in MODELS['$model']['props']))")
manifest=fd_manifest.tsv
[ -f "$manifest" ] || printf "job\tparameter\tsign\trel_step\tvalue\tstatus\tseconds\n" > "$manifest"
for h in $steps; do
  for idx in "${!names[@]}"; do
    for sign in p m; do
      job="claude_${model}_${names[$idx]}_${sign}${h/./p}"
      if grep -q "^$job	" "$manifest" 2>/dev/null; then continue; fi
      props=$(python3 - "$idx" "$sign" "$h" "${nominal[@]}" <<'EOF'
import sys
idx, sign, h = int(sys.argv[1]), sys.argv[2], float(sys.argv[3])
p = [float(x) for x in sys.argv[4:]]
p[idx] = p[idx] * (1 + h if sign == 'p' else 1 - h)
print(",".join(repr(x) for x in p))
EOF
)
      value=$(echo "$props" | cut -d, -f$((idx + 1)))
      python3 "$here/gen_cantilever.py" "$model" --props "$props" --out "$job.inp" > /dev/null
      t0=$(date +%s)
      abaqus job="$job" input="$job.inp" user="$umat" double=both cpus=4 interactive > "$job.log" 2>&1
      if grep -q "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" "$job.sta" 2>/dev/null; then
        abaqus python "$here/export_odb.py" -- "$job.odb" "${job}_fields.npz" >> "$job.log" 2>&1 && status=completed || status=export_failed
      else
        status=analysis_failed
      fi
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$job" "${names[$idx]}" "$sign" "$h" "$value" "$status" "$(( $(date +%s) - t0 ))" >> "$manifest"
      rm -f "$job.odb" "$job.com" "$job.prt" "$job.odb_f" "$job.023" "$job.mdl" "$job.stt" "$job.res" "$job.sim"
    done
  done
done
echo "done $model: $(grep -c completed "$manifest") completed of $(( $(wc -l < "$manifest") - 1 ))"
