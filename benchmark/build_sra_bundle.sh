#!/bin/bash
# Build one benchmark bundle end to end under constructed/<name>/.
#
#   benchmark/build_sra_bundle.sh <name> <N accessions> <seed> [num_workers] [from_step] [to_step]
#
# Step 1 decides which accessions are in the set and writes the frozen list
# accessions/<name>.txt; steps 2-6 only ever read that list. To build from a
# hand-made list, put it at accessions/<name>.txt and start from step 2.
#
#   1 sample exactly N accessions, excluding the training lists -> accessions/<name>.txt
#   2 download logan contigs
#   3 stream raw reads and build the initial query fasta
#   4 align queries against every accession's contigs (both strands) — expensive;
#     run on a compute node for N >~ 100 (see benchmark/align.sbatch)
#   5 finalize the bundle (filters, strand column)
#   6 write the per-rate mutated query files
set -euo pipefail
cd "$(dirname "$0")/.."

NAME=$1; N=$2; SEED=$3; WORKERS=${4:-16}; FROM=${5:-1}; TO=${6:-6}
SEQSTATS=${SEQSTATS:-/pscratch/sd/r/rsynk/rawbert_data/data/logan-seqstats-contigs-v1.1.parquet}
TABLE=accessions/TableS12_SRA_Public_100studies.tsv
EXCLUDE="[accessions/train.txt,accessions/val.txt]"
FINAL_QUERIES=${FINAL_QUERIES:-500}
NUM_READS=${NUM_READS:-1000}

ACCS=accessions/$NAME.txt
WORK=constructed/$NAME/work; BUNDLE=constructed/$NAME/bundle
mkdir -p "$WORK"

step() { [ "$FROM" -le "$1" ] && [ "$1" -le "$TO" ]; }

if step 1; then
  uv run python dataset_queries/subsample_accessions.py "$TABLE" "$SEQSTATS" "$N" "$ACCS" \
    --seed="$SEED" --exclude_files="$EXCLUDE" --max_nbseq=2000000 --min_bp=1000000
fi
[ "$(wc -l < "$ACCS")" -eq "$N" ] || { echo "$ACCS has $(wc -l < "$ACCS") accessions, expected $N"; exit 1; }
if step 2; then
  uv run python benchmark/download_logan_contigs.py "$ACCS" "$WORK/logan_contigs" --strict=true
fi
if step 3; then
  uv run python benchmark/sample_raw_read_queries.py "$ACCS" "$WORK/raw_reads" "$WORK/queries.fa" \
    --num_reads="$NUM_READS" --seed="$SEED"
fi
if step 4; then
  uv run python benchmark/get_query_alignments_minimap.py "$WORK/queries.fa" "$WORK/logan_contigs" "$WORK" \
    --num_workers="$WORKERS"
fi
if step 5; then
  uv run python benchmark/finalize_query_dataset.py "$WORK/queries.fa" "$WORK/queries_alignment_stats.parquet" \
    "$ACCS" "$BUNDLE" --final_query_num="$FINAL_QUERIES" --seed="$SEED"
fi
if step 6; then
  uv run python benchmark/mutate_queries.py --bundle "$BUNDLE"
  echo "Bundle: $BUNDLE  (accession list: $ACCS)"
fi
