# Deploying This ANNEVO Fork

This repository is a deployment-oriented private fork that includes:

- upstream ANNEVO code
- local memory-stability fixes for fragmented genomes
- bundled model weights under `ANNEVO_model/`
- exported conda environment files under `env/`

## 1. Recommended server requirements

- Linux x86_64
- Conda or Miniconda
- NVIDIA GPU recommended
- NVIDIA driver compatible with the PyTorch CUDA runtime in the exported environment
- Enough disk space for the conda environment, the repository, temporary HDF5 files, and your genomes

## 2. Clone the repository

```bash
git clone <your-private-repo-url>
cd ANNEVO_deploy_private
```

## 3. Create the conda environment

You have three environment export choices in `env/`:

- `annevo.explicit.txt`: most exact clone on a compatible Linux machine
- `annevo.no-builds.yml`: more portable conda recreation
- `annevo.full.yml`: full export including exact package builds

Preferred order:

```bash
# Most reproducible on a similar Linux server
conda create -n annevo --file env/annevo.explicit.txt

# If the explicit file is not suitable on the target machine, use:
conda env create -n annevo -f env/annevo.no-builds.yml
```

Activate the environment:

```bash
conda activate annevo
```

## 4. Verify the runtime

```bash
python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())"
python -m unittest discover -s tests -p 'test_chunking_utils.py'
python -m py_compile annotation.py prediction.py decoding.py src/predict_nucleotide.py src/chunking_utils.py
```

## 5. Important runtime note

Before running ANNEVO, this fork clears `LD_LIBRARY_PATH` to avoid CUDA or cuBLAS conflicts seen on some servers:

```bash
unset LD_LIBRARY_PATH
```

The bundled `batch.sh` already does this automatically.

## 6. Prepare input genomes

By default, `batch.sh` looks for genomes under:

```bash
./fna/*_genomic.fna.masked
```

Create the input directory and place genomes there:

```bash
mkdir -p fna
```

Outputs will be written by default to:

- `./output1/`
- `./logs/`

## 7. Run a single genome

```bash
python annotation.py \
  --genome fna/Species_genomic.fna.masked \
  --model_path ./ANNEVO_model/ANNEVO_Mammalia.pt \
  --output ./output1/Species_genomic.fna.masked.gff \
  --threads 16 \
  --genome_size_threshold 104857600 \
  --batch_size 32 \
  --max_windows_per_chunk 8192 \
  --num_workers 2
```

## 8. Run the batch script

`batch.sh` is parameterized through environment variables, so you can either use defaults or override them:

```bash
bash batch.sh
```

Example with explicit overrides:

```bash
PYTHON_BIN=python \
GENOME_DIR=./fna \
OUTDIR=./output1 \
LOGDIR=./logs \
THREADS=16 \
GENOME_SIZE_THRESHOLD=104857600 \
BATCH_SIZE=32 \
MAX_WINDOWS_PER_CHUNK=8192 \
NUM_WORKERS=2 \
bash batch.sh
```

## 9. Output naming convention

This workflow expects one genome per output set:

- input: `Species_genomic.fna.masked`
- GFF: `Species_genomic.fna.masked.gff`
- peptide FASTA: `Species_genomic.pep.fa`

For `new` genomes:

- input: `Species_new_genomic.fna.masked`
- GFF: `Species_new_genomic.fna.masked.gff`
- peptide FASTA: `Species_new_genomic.pep.fa`

## 10. Local code changes included in this fork

- streamed HDF5 prediction writes in `src/predict_nucleotide.py`
- window-cap chunking helper in `src/chunking_utils.py`
- `--max_windows_per_chunk` wired into `prediction.py` and `annotation.py`
- subprocess error propagation in `annotation.py`
- portable batch execution in `batch.sh`

## 11. Environment reference files

The `env/` directory includes:

- `annevo.explicit.txt`
- `annevo.full.yml`
- `annevo.no-builds.yml`
- `annevo.from-history.yml`
- `annevo.pip-freeze.txt`
- `runtime_versions.txt`

Use `annevo.explicit.txt` first when you want the closest possible clone of the current working environment.
