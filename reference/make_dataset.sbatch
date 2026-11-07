#!/bin/bash
#SBATCH -A m5408
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -N 1
#SBATCH -t 03:00:00
#SBATCH -c 16
#SBATCH -J locale-refcrop
#SBATCH -o constructed/slurm_logs/%j.out
#SBATCH -e constructed/slurm_logs/%j.err
# Crop the downloaded reference genomes into a training parquet on a CPU node
# (the full set is ~15 Gbp of python strings; keep it off the login node).
# Submit from the repo root:
#   NAME=reference_nohm EXCLUDE="[Homo sapiens,Mus musculus]" sbatch reference/make_dataset.sbatch
#   NAME=reference      EXCLUDE="[]"                           sbatch reference/make_dataset.sbatch
# TARGET=<rows> makes overlapping windows so the parquet has exactly that many
# rows (6768336 = the original 7-genome set); unset for non-overlapping crops.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${NAME:=reference_nohm}"
: "${EXCLUDE:=[Homo sapiens,Mus musculus]}"
: "${DATA_DIR:=/pscratch/sd/r/rsynk/reference_data/refs/ncbi_dataset/data}"
: "${TARGET:=}"
OUT=constructed/$NAME/reference_genome_dataset.parquet
mkdir -p "constructed/$NAME" constructed/slurm_logs

uv run python reference/crop_references.py "$OUT" \
    --data_dir "$DATA_DIR" --exclude_organisms "$EXCLUDE" \
    ${TARGET:+--target_rows "$TARGET"}
uv run python reference/dataset_stats.py --parquet_path "$OUT"
