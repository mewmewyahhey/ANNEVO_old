#!/usr/bin/env bash
set -euo pipefail

# Copy model_prediction.h5 files out of currently stuck one-step annotation jobs.
# This script does not stop jobs and does not delete anything.

ANNEVO_DIR="${ANNEVO_DIR:-$HOME/primate_project/ANNEVO_old}"
H5_NAME="${H5_NAME:-model_prediction.h5}"

arg_value() {
  local key="$1"
  local args="$2"
  awk -v key="$key" '
    {
      for (i = 1; i <= NF; i++) {
        if ($i == key && i < NF) {
          print $(i + 1)
          exit
        }
      }
    }
  ' <<< "$args"
}

abs_h5_path() {
  local path="$1"
  if [[ "$path" = /* ]]; then
    printf '%s\n' "$path"
  else
    printf '%s/%s\n' "$ANNEVO_DIR" "${path#./}"
  fi
}

timestamp="$(date '+%Y%m%d_%H%M%S')"

mapfile -t annotation_pids < <(pgrep -u "$USER" -f "$ANNEVO_DIR/annotation.py" || true)

if [ "${#annotation_pids[@]}" -eq 0 ]; then
  printf 'No running annotation.py processes found for ANNEVO_DIR=%s\n' "$ANNEVO_DIR"
  exit 0
fi

for ann_pid in "${annotation_pids[@]}"; do
  ann_args="$(ps -ww -p "$ann_pid" -o args= || true)"
  genome="$(arg_value "--genome" "$ann_args")"
  output="$(arg_value "--output" "$ann_args")"

  if [ -z "$output" ]; then
    printf '[WARN] annotation pid=%s has no --output argument; skipped\n' "$ann_pid" >&2
    continue
  fi

  result_dir="$(dirname "$output")"
  mkdir -p "$result_dir"

  h5_src=""
  for child_pid in $(pgrep -P "$ann_pid" || true); do
    child_args="$(ps -ww -p "$child_pid" -o args= || true)"
    child_h5="$(arg_value "--model_prediction_path" "$child_args")"
    if [ -n "$child_h5" ]; then
      h5_src="$(abs_h5_path "$child_h5")"
      break
    fi
  done

  if [ -z "$h5_src" ]; then
    printf '[WARN] annotation pid=%s has no child with --model_prediction_path; skipped\n' "$ann_pid" >&2
    continue
  fi

  if [ ! -s "$h5_src" ]; then
    printf '[WARN] source H5 missing or empty: %s\n' "$h5_src" >&2
    continue
  fi

  final_h5="$result_dir/$H5_NAME"
  done_file="$final_h5.done"

  if [ -s "$final_h5" ] && [ -f "$done_file" ]; then
    printf '[SKIP] existing rescued H5 with done marker: %s\n' "$final_h5"
    continue
  fi

  if [ -e "$final_h5" ]; then
    dest_h5="$result_dir/model_prediction.rescued.${timestamp}.pid${ann_pid}.h5"
    dest_done="$dest_h5.done"
    printf '[WARN] final H5 exists without done marker; copying to alternate path: %s\n' "$dest_h5" >&2
  else
    dest_h5="$final_h5"
    dest_done="$done_file"
  fi

  copying_h5="$dest_h5.copying.${timestamp}.$$"
  printf '[COPY] pid=%s\n' "$ann_pid"
  printf '       source=%s\n' "$h5_src"
  printf '       dest=%s\n' "$dest_h5"
  cp -a "$h5_src" "$copying_h5"
  mv "$copying_h5" "$dest_h5"

  {
    printf 'annotation_pid=%s\n' "$ann_pid"
    printf 'genome=%s\n' "$genome"
    printf 'output=%s\n' "$output"
    printf 'source_h5=%s\n' "$h5_src"
    printf 'rescued_h5=%s\n' "$dest_h5"
    printf 'rescued_at=%s\n' "$(date '+%F %T')"
  } > "$dest_done"
  printf '[DONE] %s\n' "$dest_h5"
done
