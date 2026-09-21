---
pretty_name: LOCALE Benchmark - SRA-50 (v2)
license: cc-by-4.0
size_categories:
- n<1K
task_categories:
- feature-extraction
- sentence-similarity
tags:
- genomics
- dna
- bioinformatics
- vector-search
- sequence-search
- retrieval
- sra
configs:
- config_name: queries
  data_files: queries.parquet
---

# LOCALE Benchmark: SRA-50 (v2)

Given a sequencing read, find the SRA runs whose assembled contigs it aligns to. This is the 50-run set: **50 SRA runs** indexed, **500 query reads**, ground truth from minimap2 alignment of each read against every run's [Logan](https://github.com/IndexThePlanet/Logan) contigs. It is the small-scale companion of [SRA-500 v2](https://huggingface.co/datasets/rsynk/locale-benchmark-sra500-v2) and [SRA-4571](https://huggingface.co/datasets/rsynk/locale-benchmark-sra4571), built by the same pipeline with the same filters.

This release supersedes [`rsynk/locale-benchmark-sra50`](https://huggingface.co/datasets/rsynk/locale-benchmark-sra50) (47 runs). The runs are a fresh seeded draw (seed 50) and share no accessions with the old set. The old set's queries were filtered so that the run a read came from was *absent* from its relevant set, an inverted condition; here every query's source run is in its relevant set, alignment ran on both strands, and a `strand` column was added. Embeddings for this set are not yet published.

| | |
|---|---|
| Runs indexed | 50 |
| Queries | 500 (from 43 of the runs) |
| Relevant runs per query | mean 2.51, median 1, max 10 |
| Query length | 51 to 170 bp, mean 128 |
| Queries on the `-` strand | 276 (55%) |
| Filter accounting | 50000 sampled reads; 35114 dropped for no hit to their own run, 1220 for hitting more than 10 runs, 7777 for overlapping another query; 5889 eligible, 500 kept |
| Seed | 50 |
| Excluded from the draw | the LOCALE training and validation runs (521) |

## How the ground truth was built

1. **Accession set.** Drawn from the SRA runs of the [Metagraph 100-studies table](https://github.com/ratschlab/metagraph) proportionally by organism, restricted to runs with at most 2M Logan contigs and at least 1 Mbp of assembled sequence, certified to exist in Logan at draw time, and with every run in the LOCALE training and validation sets removed from the pool before the draw. The draw is seeded, so `accs.txt` is reproducible from the seed.
2. **Queries.** The first 1000 spots of every run were streamed from SRA with `fastq-dump` (mate 1 for paired runs). Stored query sequences are the raw reads exactly as sequenced; nothing is trimmed or reverse-complemented.
3. **Alignment.** Every read was aligned with minimap2 (`-a --eqx`, both strands) against the Logan contigs of every run in the set, keeping the best alignment per (read, run) by exact-match count. An alignment counts when identity (exact matches / read length) exceeds 0.9.
4. **Filters**, per query, in this order: the run the read was sampled from must be in its relevant set; the read must not hit more than 20% of the runs in the set; the read must not share a contig region with another query (strict overlap filter with a 1024 bp chunk margin). Survivors are shuffled with the seed and the first N kept.
5. **Mutated copies.** The same queries with sequencing noise injected by [mutation-simulator](https://github.com/mkpython3/Mutation-Simulator) at SNP rate r with insertion and deletion rates r/10 each, r = 0.05 and 0.10. Every row keeps its `query_id` and ground-truth columns, so any method can be evaluated on identical mutated reads. The simulator's fasta and VCF per rate are under `mutations/`.

## Files

| File | Contents |
|---|---|
| `accs.txt` | the indexed SRA run accessions, one per line |
| `queries.parquet` | one row per query with its ground truth (schema below) |
| `queries_mut0.00.parquet` | `queries.parquet` plus a `mutation_rate` column (clean) |
| `queries_mut0.05.parquet`, `queries_mut0.10.parquet` | same rows in the same order with `query_sequence` mutated |
| `queries.fa`, `queries.fa.fai` | the clean queries as fasta (mutation-simulator's input) and its index |
| `mutations/` | mutation-simulator fasta and VCF for each rate |

## Schema of `queries.parquet`

| Column | Type | Meaning |
|---|---|---|
| `accession` | str | run the read was sampled from |
| `query_id` | str | `<accession>.<spot>` from fastq-dump |
| `query_sequence` | str | the raw read |
| `contig_accession` | list[str] | **relevant set**: every run with an alignment above the identity threshold; always contains `accession` |
| `contig_id` | list[str] | best-aligned contig in each relevant run |
| `identity` | list[f64] | exact matches / read length of that alignment |
| `ratio` | list[f64] | contig length / read length |
| `contig_len` | list[i64] | contig length |
| `aln_interval_contig` | list[struct] | 0-based `[start, end)` of the alignment on the contig |
| `strand` | str | `+`/`-`: strand of the best hit to the source run |

The six list columns are parallel: element *i* of each describes the same alignment.
## License

The data is licensed under Creative Commons `cc-by-4.0`; the code that built it is MIT ([locale-data](https://github.com/ryansynk/locale-data)). The underlying sequence data comes from the NCBI Sequence Read Archive, which places [no restrictions on redistribution](https://www.ncbi.nlm.nih.gov/home/about/policies/); assembled contigs come from [Logan](https://github.com/IndexThePlanet/Logan).

## Citation

For questions, reach out to `ryansynk@umd.edu`. If you find this work useful, please cite:

```
@article {Synk2026.05.12.724581,
	author = {Synk, Ryan P. and Pandey, Prashant and Sahinalp, S. Cenk and Duraiswami, Ramani},
	title = {LOCALE: Local-Alignment Embeddings for Noise-Robust DNA Search at SRA Scale},
	elocation-id = {2026.05.12.724581},
	year = {2026},
	doi = {10.64898/2026.05.12.724581},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2026/05/14/2026.05.12.724581},
	eprint = {https://www.biorxiv.org/content/early/2026/05/14/2026.05.12.724581.full.pdf},
	journal = {bioRxiv}
}
```
