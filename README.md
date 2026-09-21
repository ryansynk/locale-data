# locale-data

Scripts for generating the training and benchmark datasets for the locale encoder model.

```
dataset_queries/   SQL against SRA metadata (BigQuery) + subsample script → candidate accession lists
accessions/        frozen accession lists (the exact lists behind the published HF datasets)
benchmark/         candidate list → benchmark bundle (queries.parquet, queries_mut<rate>.parquet, queries.fa, accs.txt)
training/          candidate list → training parquet (all_contigs.parquet)
constructed/       work dirs and bundles, `constructed/<set>/{work,bundle}` (gitignored)
```

Both pipelines start from an **accession list** (one accession per line). Deciding what is in the list (organism proportions, size filters, existence in logan, exclusion of training accessions) happens entirely in the selection step; the build steps only read the list. `benchmark/download_logan_contigs.py` also emits the list it actually downloaded (`downloaded_accs.txt`) for training draws where a few candidates may be missing from logan.

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
    500 accessions/sra500.txt --seed=500 \
    --exclude_files="[accessions/train.txt,accessions/val.txt]" \
    --max_nbseq=2000000 --min_bp=1000000
```

Keep train, val, and benchmark lists disjoint (in SQL, or via `--exclude_files=[...]`, which removes those accessions from the pool before the draw). The draw is deterministic in `--seed`, and each candidate is checked against logan-pub as it is drawn so the list has exactly N accessions (`--verify_logan=false` to skip for large training draws). The lists in `accessions/` are frozen: they are the exact draws behind the published datasets.

## Benchmark pipeline

Produces a bundle in the layout of the Hugging Face benchmark datasets, ready to upload as-is. `benchmark/build_sra_bundle.sh <set> <N> <seed>` runs the draw above and the steps below under `constructed/<set>/` (it numbers the draw as step 1, so `from_step`/`to_step` are offset by one from the list below); `benchmark/align.sbatch` runs the alignment onwards on a compute node.

```bash
ACCS=accessions/<set>.txt
WORK=constructed/<set>/work

# 1. Download logan contigs. --strict: the list was checked at draw time,
#    so a missing accession is a real error.
uv run python benchmark/download_logan_contigs.py $ACCS $WORK/logan_contigs --strict=true

# 2. Stream the first N raw reads per accession from SRA (no full download)
#    and build the initial query fasta.
uv run python benchmark/sample_raw_read_queries.py $ACCS $WORK/raw_reads $WORK/queries.fa --num_reads=1000 --seed=$SEED

# 3. Align every query against every accession's contigs, on both strands.
#    Writes $WORK/alignments/<acc>.sam and $WORK/queries_alignment_stats.parquet
#    (parquet is named after the query fasta's stem). Expensive — use a compute node.
uv run python benchmark/get_query_alignments_minimap.py $WORK/queries.fa $WORK/logan_contigs $WORK --num_workers=32

# 4. Filter alignments (identity > 0.9; then per query: source accession must
#    be in the relevant set, max-match, overlap) and sample the final queries
#    (seeded) into a self-contained bundle.
uv run python benchmark/finalize_query_dataset.py \
    $WORK/queries.fa $WORK/queries_alignment_stats.parquet $ACCS $WORK/bundle \
    --final_query_num=500 --seed=$SEED

# 5. Write one pre-mutated query file per mutation rate (0.00 / 0.05 / 0.10)
#    with mutation-simulator: SNP rate = rate, insertion and deletion rates =
#    rate/10 each. The benchmark never mutates on the fly; its --mutation_rate
#    picks one of these files.
uv run python benchmark/mutate_queries.py --bundle $WORK/bundle
```

`$WORK/bundle/` then contains:

- `queries.parquet` — the clean queries with their alignment ground truth; what `print_results.py` and `make_table3.py` read. `query_sequence` is the raw read as sequenced (never reverse-complemented); `strand` (`+`/`-`) is the strand of the best hit to the source accession and is the only column added to the published schema.
- `queries_mut0.00.parquet`, `queries_mut0.05.parquet`, `queries_mut0.10.parquet` — `queries.parquet` with `query_sequence` replaced by the mutated read and a `mutation_rate` column added. Same rows in the same order, so the benchmark's seeded subsample picks the same queries at every rate. `run_benchmark.py --mutation_rate <rate>` loads `queries_mut<rate>.parquet`.
- `queries.fa` — the clean queries as fasta (mutation-simulator's input).
- `mutations/` — mutation-simulator's fasta and VCF per rate. mutation-simulator has no seed flag, so these are the provenance record: rerunning step 5 draws different mutations, and the parquets are what the benchmark reads.
- `accs.txt` — a copy of the true list, with an assertion that the queries only reference accessions in it.

Pass `--rates=[0.0,0.02]` or `--indel_fraction=0.0` to step 5 to change the rates; any rate the benchmark is run at must have its file in the bundle.

### SRA-viral (HBV) bundle

`benchmark/build_viral_bundle.py` builds the cross-genotype retrieval bundle: the index set is an accession list (the sra50 distractors) plus the genotype-B runs (`accessions/hbv_genotype_B.txt`), queries are raw reads from the genotype-D runs (`accessions/hbv_genotype_D.txt`), and every query's relevant set is the five genotype-B accessions (Rq = 5) — relevance comes from SRA genotype metadata, not minimap2. Raw-read sampling and mutation are the same scripts as above; there is no alignment step.

```bash
uv run python benchmark/build_viral_bundle.py \
    --index_accs accessions/sra50.txt \
    --relevant_accs accessions/hbv_genotype_B.txt \
    --query_accs accessions/hbv_genotype_D.txt \
    --work_dir constructed/sra55viral/work --output_path constructed/sra55viral/bundle \
    --seed=55
```

The bundle has the standard layout and schema (`identity`, `ratio`, `contig_len`, `aln_interval_contig` and `strand` are null; `contig_id` holds the relevant accession), so `run_benchmark.py` runs on it unchanged. `accessions/sra55viral.txt` is its frozen index list (sra50 + genotype B).

## Training pipeline

The locale trainer (`data_type: contig`) consumes one parquet per split with a `sequence` column; cropping, mutation, and length filtering all happen at train time in the Batcher/Augmenter. Building a dataset is therefore just: download contigs, chunk, shuffle, write parquet.

```bash
# 1. Download (same script as the benchmark pipeline); missing accessions are
#    reported, not fatal, and downloaded_accs.txt is the list actually built
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
