#!/usr/bin/env bash
#SBATCH --job-name=annevo-a100-01
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=annevo-a100-01-%j.out
#SBATCH --error=annevo-a100-01-%j.err

set -uo pipefail
unset LD_LIBRARY_PATH

BATCH_ID=1
TOTAL_BATCHES=8

activate_conda_env() {
  set +u
  conda activate "$1"
  set -u
}

ANNEVO_DIR="${ANNEVO_DIR:-$HOME/primate_project/ANNEVO_old}"
DATA_ROOT="${DATA_ROOT:-$HOME/primate_project/annotation_process/data}"
ANNEVO_ENV="${ANNEVO_ENV:-ANNEVO}"
CONDA_BASE="${CONDA_BASE:-$HOME/miniconda3}"
LIST_PATH="${LIST_PATH:-$ANNEVO_DIR/tba_list.txt}"
MODEL_PATH="${MODEL_PATH:-$ANNEVO_DIR/ANNEVO_model/ANNEVO_Mammalia.pt}"
LINEAGE="${LINEAGE:-Mammalia}"
GENOME_SIZE_THRESHOLD="${GENOME_SIZE_THRESHOLD:-104857600}"
BATCH_SIZE="${BATCH_SIZE:-16}"
MAX_WINDOWS_PER_CHUNK="${MAX_WINDOWS_PER_CHUNK:-8192}"
NUM_WORKERS="${NUM_WORKERS:-2}"
OVERLAP_PRED="${OVERLAP_PRED:-0}"
SKIP_DONE="${SKIP_DONE:-1}"
H5_NAME="${H5_NAME:-model_prediction.h5}"
TMP_ROOT="${TMP_ROOT:-$ANNEVO_DIR/tmp}"

cd "$ANNEVO_DIR" || exit 1

if [ ! -f "$LIST_PATH" ]; then
  printf 'ERROR: species list not found: %s\n' "$LIST_PATH" >&2
  exit 1
fi

if [ ! -f "$MODEL_PATH" ]; then
  printf 'ERROR: ANNEVO model not found: %s\n' "$MODEL_PATH" >&2
  exit 1
fi

if command -v module >/dev/null 2>&1; then
  module load cuda/12.1.1 || true
fi

source "$CONDA_BASE/etc/profile.d/conda.sh"
activate_conda_env "$ANNEVO_ENV"

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="$NUM_WORKERS"
export TOKENIZERS_PARALLELISM=false

printf 'job_id=%s\n' "${SLURM_JOB_ID:-unknown}"
printf 'node=%s\n' "${SLURMD_NODENAME:-unknown}"
printf 'batch=%s/%s\n' "$BATCH_ID" "$TOTAL_BATCHES"
printf 'annevo_dir=%s\n' "$ANNEVO_DIR"
printf 'data_root=%s\n' "$DATA_ROOT"
printf 'model=%s\n' "$MODEL_PATH"
printf 'batch_size=%s genome_size_threshold=%s num_workers=%s overlap_pred=%s\n' "$BATCH_SIZE" "$GENOME_SIZE_THRESHOLD" "$NUM_WORKERS" "$OVERLAP_PRED"

if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi
fi

python -u -c "import torch; print('torch=', torch.__version__, 'cuda=', torch.version.cuda, 'cuda_available=', torch.cuda.is_available(), flush=True)"

mapfile -t ALL_SPECIES < <(awk 'NF {print $1}' "$LIST_PATH")
species_total="${#ALL_SPECIES[@]}"
base_count=$((species_total / TOTAL_BATCHES))
remainder=$((species_total % TOTAL_BATCHES))
if [ "$BATCH_ID" -le "$remainder" ]; then
  start_index=$(((BATCH_ID - 1) * (base_count + 1)))
  batch_count=$((base_count + 1))
else
  start_index=$((remainder * (base_count + 1) + (BATCH_ID - remainder - 1) * base_count))
  batch_count=$base_count
fi
SPECIES_LIST=("${ALL_SPECIES[@]:start_index:batch_count}")

printf 'species_total=%s start_line=%s end_line=%s\n' "$species_total" "$((start_index + 1))" "$((start_index + batch_count))"
printf 'species_count=%s\n' "${#SPECIES_LIST[@]}"

success_count=0
skip_count=0
fail_count=0

for species in "${SPECIES_LIST[@]}"; do
  genome="$DATA_ROOT/$species/genome/${species}_genomic.fna.masked"
  result_dir="$DATA_ROOT/$species/results/09_annevo"
  output_h5="$result_dir/$H5_NAME"
  done_file="$output_h5.done"

  printf '\n[%s] START species=%s\n' "$(date '+%F %T')" "$species"
  printf 'genome=%s\n' "$genome"
  printf 'output=%s\n' "$output_h5"

  mkdir -p "$result_dir"

  done_markers=("$result_dir"/*.done)
  h5_files=("$result_dir"/*.h5)
  if [ "$SKIP_DONE" = "1" ] && { [ -e "${done_markers[0]}" ] || [ -e "${h5_files[0]}" ]; }; then
    printf '[%s] SKIP existing done/H5 result: %s\n' "$(date '+%F %T')" "$species"
    skip_count=$((skip_count + 1))
    continue
  fi

  if [ ! -s "$genome" ]; then
    printf '[%s] ERROR missing genome: %s\n' "$(date '+%F %T')" "$genome" >&2
    fail_count=$((fail_count + 1))
    continue
  fi

  run_dir="$TMP_ROOT/annevo_prediction/$species/${SLURM_JOB_ID:-manual}.$(date '+%Y%m%d_%H%M%S')"
  mkdir -p "$run_dir"
  run_h5="$run_dir/$H5_NAME"

  cmd=(
    python -u "$ANNEVO_DIR/prediction.py"
    --genome "$genome"
    --model_path "$MODEL_PATH"
    --model_prediction_path "$run_h5"
    --genome_size_threshold "$GENOME_SIZE_THRESHOLD"
    --batch_size "$BATCH_SIZE"
    --max_windows_per_chunk "$MAX_WINDOWS_PER_CHUNK"
    --num_workers "$NUM_WORKERS"
  )

  printf '[%s] CMD:' "$(date '+%F %T')"
  printf ' %q' "${cmd[@]}"
  printf '\n'

  if "${cmd[@]}"; then
    if [ ! -s "$run_h5" ]; then
      printf '[%s] FAIL species=%s no non-empty H5 produced: %s\n' "$(date '+%F %T')" "$species" "$run_h5" >&2
      fail_count=$((fail_count + 1))
      continue
    fi
    cp -a "$run_h5" "$output_h5"
    printf 'species=%s\nfinished_at=%s\noutput=%s\nrun_h5=%s\n' "$species" "$(date '+%F %T')" "$output_h5" "$run_h5" > "$done_file"
    printf '[%s] DONE species=%s\n' "$(date '+%F %T')" "$species"
    success_count=$((success_count + 1))
  else
    status=$?
    printf '[%s] FAIL species=%s exit_code=%s\n' "$(date '+%F %T')" "$species" "$status" >&2
    fail_count=$((fail_count + 1))
  fi
done

printf '\n[%s] SUMMARY batch=%s success=%s skipped=%s failed=%s\n' "$(date '+%F %T')" "$BATCH_ID" "$success_count" "$skip_count" "$fail_count"

if [ "$fail_count" -gt 0 ]; then
  exit 1
fi
