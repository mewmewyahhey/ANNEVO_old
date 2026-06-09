#!/usr/bin/env bash
unset LD_LIBRARY_PATH
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MODEL="${MODEL:-$SCRIPT_DIR/ANNEVO_model/ANNEVO_Mammalia.pt}"
PYTHON_BIN="${PYTHON_BIN:-python}"
THREADS="${THREADS:-16}"
GENOME_DIR="${GENOME_DIR:-$SCRIPT_DIR/fna}"
OUTDIR="${OUTDIR:-$SCRIPT_DIR/output1}"
LOGDIR="${LOGDIR:-$SCRIPT_DIR/logs}"
GENOME_SIZE_THRESHOLD="${GENOME_SIZE_THRESHOLD:-104857600}"
BATCH_SIZE="${BATCH_SIZE:-32}"
MAX_WINDOWS_PER_CHUNK="${MAX_WINDOWS_PER_CHUNK:-8192}"
NUM_WORKERS="${NUM_WORKERS:-2}"
INPUT_GLOB="${INPUT_GLOB:-*_genomic.fna.masked}"

mkdir -p "$OUTDIR" "$LOGDIR"

if [ ! -f "$MODEL" ]; then
  echo "Model file not found: $MODEL"
  exit 1
fi

shopt -s nullglob
files=("$GENOME_DIR"/$INPUT_GLOB)
shopt -u nullglob

if [ ${#files[@]} -eq 0 ]; then
  echo "No genome files found under $GENOME_DIR (pattern: $INPUT_GLOB)"
  exit 1
fi

for genome in "${files[@]}"; do
  base="$(basename "$genome")"

  out_gff="${OUTDIR}/${base}.gff"
  log="${LOGDIR}/${base}.log"

  if [ -s "$out_gff" ]; then
    echo "[$(date '+%F %T')] SKIP   $base (exists: $out_gff)"
    continue
  fi

  echo "[$(date '+%F %T')] START  $genome -> $out_gff"
  if "$PYTHON_BIN" annotation.py \
      --genome "$genome" \
      --model_path "$MODEL" \
      --output "$out_gff" \
      --threads "$THREADS" \
      --genome_size_threshold "$GENOME_SIZE_THRESHOLD" \
      --batch_size "$BATCH_SIZE" \
      --max_windows_per_chunk "$MAX_WINDOWS_PER_CHUNK" \
      --num_workers "$NUM_WORKERS" \
      >"$log" 2>&1
  then
      echo "[$(date '+%F %T')] DONE   $base"
  else
      rc=$?
      echo "[$(date '+%F %T')] FAIL   $base (exit=$rc). See $log"
      continue
  fi
done
