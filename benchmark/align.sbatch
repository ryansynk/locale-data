#!/bin/bash
#SBATCH -A m5408
#SBATCH -C cpu
#SBATCH -q premium
#SBATCH -N 1
#SBATCH -t 08:00:00
#SBATCH -c 256
#SBATCH -J locale-align
#SBATCH -o constructed/slurm_logs/%j.out
#SBATCH -e constructed/slurm_logs/%j.err
# Steps 4-6 of build_sra_bundle.sh on one CPU node (the alignment is the
# expensive step). Submit from the repo root after steps 1-3:
#   sbatch --export=ALL,NAME=sra500,N=500,SEED=500 benchmark/align.sbatch
cd "$SLURM_SUBMIT_DIR"
benchmark/build_sra_bundle.sh "$NAME" "$N" "$SEED" 128 4
