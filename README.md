# locale-data

Scripts for generating the training and benchmark datasets for the locale encoder model.

```
dataset_queries/   SQL against SRA metadata (BigQuery) + subsample script → candidate accession lists
accessions/        frozen accession lists (the exact lists behind the published HF datasets)
benchmark/         candidate list → benchmark bundle (queries.parquet, queries_mut<rate>.parquet, queries.fa, accs.txt)
training/          candidate list → training parquet (all_contigs.parquet)
```

Both pipelines start from a **candidate accession list** (one accession per line) and share the same first step: `benchmark/download_logan_contigs.py` certifies the candidates against what actually exists in logan and emits the **true list** (`downloaded_accs.txt`), which is what flows through everything downstream.

## Dependencies

On PATH:
- [aws-cli](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html)
- [sra-toolkit](https://github.com/ncbi/sra-tools) (`fastq-dump`)
- [minimap2](https://github.com/lh3/minimap2)

Python: `uv sync` (or just prefix commands with `uv run`).

## Selecting accessions

Candidate lists come from `dataset_queries/paired_illumina_ont.sql` run against the SRA metadata on BigQuery (`nih-sra-datastore.sra.metadata`), or by sampling from the metagraph 100-studies table:

```bash
uv run python dataset_queries/subsample_accessions.py \
    accessions/TableS12_SRA_Public_100studies.tsv \
    /path/to/logan-seqstats-contigs-v1.1.parquet \
    500 candidates.txt \
    --max_nbseq=2000000 --min_bp=1000000
```

Keep train, val, and benchmark lists disjoint (in SQL, or via `--exclude_files=[...]`). The lists in `accessions/` are frozen: they are the exact draws behind the published datasets.

## Benchmark pipeline

Produces a bundle in the layout of the Hugging Face benchmark datasets, ready to upload as-is.

```bash
WORK=/path/to/workdir

# 1. Download logan contigs; emits the true accession list.
#    Missing accessions are reported, not fatal (--strict=true to change that).
uv run python benchmark/download_logan_contigs.py candidates.txt $WORK/logan_contigs
ACCS=$WORK/logan_contigs/downloaded_accs.txt

# 2. Stream the first N raw reads per accession from SRA (no full download)
#    and build the initial query fasta.
uv run python benchmark/sample_raw_read_queries.py $ACCS $WORK/raw_reads $WORK/queries.fa --num_reads=1000

# 3. Align every query against every accession's contigs.
#    Writes $WORK/alignments/<acc>.sam and $WORK/queries_alignment_stats.parquet
#    (parquet is named after the query fasta's stem). Expensive — use a compute node.
uv run python benchmark/get_query_alignments_minimap.py $WORK/queries.fa $WORK/logan_contigs $WORK --num_workers=32

# 4. Filter alignments (identity > 0.9, overlap and max-match filters) and
#    sample the final queries into a self-contained bundle.
uv run python benchmark/finalize_query_dataset.py \
    $WORK/queries.fa $WORK/queries_alignment_stats.parquet $ACCS $WORK/bundle \
    --final_query_num=500

# 5. Write one pre-mutated query file per mutation rate (0.00 / 0.05 / 0.10)
#    with mutation-simulator: SNP rate = rate, insertion and deletion rates =
#    rate/10 each. The benchmark never mutates on the fly; its --mutation_rate
#    picks one of these files.
uv run python benchmark/mutate_queries.py --bundle $WORK/bundle
```

`$WORK/bundle/` then contains:

- `queries.parquet` — the clean queries with their alignment ground truth; what `print_results.py` and `make_table3.py` read.
- `queries_mut0.00.parquet`, `queries_mut0.05.parquet`, `queries_mut0.10.parquet` — `queries.parquet` with `query_sequence` replaced by the mutated read and a `mutation_rate` column added. Same rows in the same order, so the benchmark's seeded subsample picks the same queries at every rate. `run_benchmark.py --mutation_rate <rate>` loads `queries_mut<rate>.parquet`.
- `queries.fa` — the clean queries as fasta (mutation-simulator's input).
- `mutations/` — mutation-simulator's fasta and VCF per rate. mutation-simulator has no seed flag, so these are the provenance record: rerunning step 5 draws different mutations, and the parquets are what the benchmark reads.
- `accs.txt` — a copy of the true list, with an assertion that the queries only reference accessions in it.

Pass `--rates=[0.0,0.02]` or `--indel_fraction=0.0` to step 5 to change the rates; any rate the benchmark is run at must have its file in the bundle.

## Training pipeline

The locale trainer (`data_type: contig`) consumes one parquet per split with a `sequence` column; cropping, mutation, and length filtering all happen at train time in the Batcher/Augmenter. Building a dataset is therefore just: download contigs, chunk, shuffle, write parquet.

```bash
# 1. Certify against logan and download (same script as the benchmark pipeline)
uv run python benchmark/download_logan_contigs.py train_candidates.txt $WORK/train/logan_contigs
uv run python benchmark/download_logan_contigs.py val_candidates.txt $WORK/val/logan_contigs

# 2. Chunk (>1024bp split into 1024bp pieces), filter (<200bp dropped),
#    shuffle, and write the training parquet
uv run python training/collect_contigs.py $WORK/train/logan_contigs $WORK/train/all_contigs.parquet
uv run python training/collect_contigs.py $WORK/val/logan_contigs $WORK/val/all_contigs.parquet
```

Point the locale training config at the outputs:

```yaml
data_type: contig
dataset_path: $WORK/train/all_contigs.parquet
val_dataset_path: $WORK/val/all_contigs.parquet
```

`max_contig_len` (default 1024) should match `augment_config.max_seq_len`, and `min_contig_len` (default 200) must be >= `augment_config.min_seq_len` (150) or the Batcher will just refilter at load.
