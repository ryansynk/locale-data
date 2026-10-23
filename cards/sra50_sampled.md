---
pretty_name: LOCALE Benchmark - SRA-50 sampled
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

# LOCALE Benchmark: SRA-50 sampled

A reduced-index copy of [SRA-50 v2](https://huggingface.co/datasets/rsynk/locale-benchmark-sra50-v2) for cheap, repeated evaluation, such as tracking a model during training. The **50 runs, 500 queries, and ground truth are identical** to SRA-50 v2; only the contigs to index are subsampled, so the whole set embeds in minutes instead of embedding all 8.7M contigs. Scores are **not comparable** to SRA-50 v2: dropping most contigs raises the ratio of relevant to irrelevant vectors, so absolute numbers come out higher. Use it to compare models against each other, and SRA-50 v2 or [SRA-500 v2](https://huggingface.co/datasets/rsynk/locale-benchmark-sra500-v2) to report.

## Which contigs are kept

Every contig of every run is kept if it is either of:

1. **A relevant contig.** The best-aligned contig of every (query, run) pair in the ground truth, so no query loses a positive and `queries.parquet` is unchanged.
2. **A hit contig.** Any contig that held a top-100 window for any query in an exact search of SRA-50 v2 with LOCALE (checkpoint `8vqiabk9`, 256 bp windows, stride 106) at mutation rates 0.05 and 0.10. Most of these are wrong hits, the near-homologs a model has to rank below the true run, and a uniform subsample would discard nearly all of them. They are hard for that one model; every model evaluated here sees the same set.

Everything else is kept with probability `subsample_fraction` by a seeded coin flip per contig, so each run keeps about that share of its contigs.

| | |
|---|---|
| Runs indexed | 50 (all of SRA-50 v2) |
| Queries | 500 (all of SRA-50 v2) |
| Contigs in SRA-50 v2 | 8,712,391 |
| Relevant contigs kept | 1,250 |
| Hit contigs kept | 84,308 (1,089 of them also relevant) |
| Random fraction | 0.1 |
| Seed | 0 |
| Contigs kept | 946,972 (10.9% of SRA-50 v2; 84,469 by rule, the rest by coin flip) |

Built by `scripts/build_subsample_dataset.py` in [locale-data](https://github.com/ryansynk/locale-data); the hit contigs came from `benchmark/scripts/hit_contigs.py` in [locale](https://github.com/ryansynk/locale) run on the SRA-50 v2 exact-search hits.

## Files

The layout is SRA-50 v2's without its `mutations/` provenance directory. Every file except the contigs is a verbatim copy.

| File | Contents |
|---|---|
| `accs.txt` | the indexed SRA run accessions, one per line |
| `logan_accessions/<run>/<run>.contigs.fa` | the kept contigs of each run, in their original order |
| `queries.parquet` | one row per query with its ground truth (schema below) |
| `queries_mut0.00.parquet` | `queries.parquet` plus a `mutation_rate` column (clean) |
| `queries_mut0.05.parquet`, `queries_mut0.10.parquet` | same rows in the same order with `query_sequence` mutated |
| `queries.fa`, `queries.fa.fai` | the clean queries as fasta (mutation-simulator's input) and its index |

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

The six list columns are parallel: element *i* of each describes the same alignment. See the SRA-50 v2 card for how the queries and ground truth were built.

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
