#!/usr/bin/env bash
set -euo pipefail

script="${1:-batch.sh}"

if [ ! -f "$script" ]; then
  printf 'ERROR: job script not found: %s\n' "$script" >&2
  exit 1
fi

backup="${script}.before_h5_prediction.$(date '+%Y%m%d_%H%M%S')"
cp -a "$script" "$backup"

# Normalize CRLF first; otherwise bash may fail on scripts copied from Windows.
sed -i 's/\r$//' "$script"

grep -q '^MAX_WINDOWS_PER_CHUNK=' "$script" || \
  sed -i '/^BATCH_SIZE=/a MAX_WINDOWS_PER_CHUNK="${MAX_WINDOWS_PER_CHUNK:-8192}"' "$script"

grep -q '^H5_NAME=' "$script" || \
  sed -i '/^SKIP_DONE=/a H5_NAME="${H5_NAME:-model_prediction.h5}"' "$script"

grep -q '^TMP_ROOT=' "$script" || \
  sed -i '/^H5_NAME=/a TMP_ROOT="${TMP_ROOT:-$ANNEVO_DIR/tmp}"' "$script"

sed -i -f <(cat <<'SED'
s#^  output_gff="\$result_dir/\${species}_genomic\.fna\.masked\.gff"$#  output_h5="$result_dir/$H5_NAME"#
s#^  done_file="\$output_gff\.done"$#  done_file="$output_h5.done"#
s#^  printf 'output=%s\\n' "\$output_gff"$#  printf 'output=%s\\n' "$output_h5"#

/^  if \[ "\$SKIP_DONE" = "1" \] && \[ -s "\$output_gff" \] && \[ -f "\$done_file" \]; then$/,/^    printf '\[%s\] SKIP existing completed result: %s\\n' "\$(date '+%F %T')" "\$species"$/c\
  done_markers=("$result_dir"/*.done); h5_files=("$result_dir"/*.h5); if [ "$SKIP_DONE" = "1" ] && { [ -e "${done_markers[0]}" ] || [ -e "${h5_files[0]}" ]; }; then printf '[%s] SKIP existing done/H5 result: %s\\n' "$(date '+%F %T')" "$species"

/^  if \[ -e "\$output_gff" \] && \[ ! -f "\$done_file" \]; then$/,/^  fi$/c\
  run_dir="$TMP_ROOT/annevo_prediction/$species/${SLURM_JOB_ID:-manual}.$(date '+%Y%m%d_%H%M%S')"; mkdir -p "$run_dir"; run_h5="$run_dir/$H5_NAME"

s#^    python -u "\$ANNEVO_DIR/annotation\.py"$#    python -u "$ANNEVO_DIR/prediction.py"#
s#^    --output "\$output_gff"$#    --model_prediction_path "$run_h5"#

/^    --threads 16$/,/^    --num_workers 2$/c\
    --max_windows_per_chunk "$MAX_WINDOWS_PER_CHUNK" --num_workers "$NUM_WORKERS"

/^  if \[ "\$OVERLAP_PRED" = "1" \]; then$/,/^  fi$/d

/^    printf 'species=%s\\nfinished_at=%s\\noutput=%s\\n' "\$species" "\$(date '+%F %T')" "\$output_gff" > "\$done_file"$/c\
    if [ ! -s "$run_h5" ]; then printf '[%s] FAIL species=%s no non-empty H5 produced: %s\\n' "$(date '+%F %T')" "$species" "$run_h5" >&2; fail_count=$((fail_count + 1)); continue; fi; cp -a "$run_h5" "$output_h5"; printf 'species=%s\\nfinished_at=%s\\noutput=%s\\nrun_h5=%s\\n' "$species" "$(date '+%F %T')" "$output_h5" "$run_h5" > "$done_file"
SED
) "$script"

printf 'converted=%s\nbackup=%s\n' "$script" "$backup"
