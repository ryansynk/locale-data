---
pretty_name: LOCALE Benchmark - SRA-55 cross-genotype viral retrieval
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
- hepatitis-b-virus
configs:
- config_name: queries
  data_files: queries.parquet
---

# LOCALE Benchmark: SRA-55 cross-genotype viral retrieval

A retrieval benchmark built on **real biological divergence** between viral strains rather than noise injected into reads. Query with sequencing reads from one Hepatitis B virus genotype; retrieve the SRA runs of a *different* genotype. The two genotypes are about 11% divergent, enough to break exact k-mer overlap while remaining genuinely homologous.

This release supersedes [`rsynk/locale-benchmark-sra52-viral`](https://huggingface.co/datasets/rsynk/locale-benchmark-sra52-viral). The relevant and query runs are unchanged; the distractors are now the 50 runs of [SRA-50 v2](https://huggingface.co/datasets/rsynk/locale-benchmark-sra50-v2) instead of the 47 of the old SRA-50, and the schema matches the alignment benchmarks exactly.

| | |
|---|---|
| Runs indexed | **55** (50 SRA-50 v2 distractors + 5 relevant) |
| Queries | **500** reads |
| Query genotype | **D** (10 runs, BioProject PRJNA659786) |
| Relevant genotype | **B** (5 runs, PRJNA556730 and PRJNA564199) |
| Relevance density | 5 / 55 = 9.1% |
| Measured divergence between the genotypes | 0.890 identity, 0.983 coverage |
| Query length | 35 to 151 bp, mean 118 |
| Seed | 55 |

## Ground truth is metadata, not alignment

`contig_accession` holds the same five genotype-B runs for every query: relevance comes from the SRA genotype annotation, so no aligner decides the correct answer. The alignment-derived columns exist so the file has the same schema as the SRA-50/500/4571 benchmarks and the same code runs on it, but `identity`, `ratio`, `contig_len`, `aln_interval_contig` and `strand` are null, and `contig_id` holds the relevant run accession itself. The query runs are not indexed.

## How it was built

1. **Index set.** The 50 runs of SRA-50 v2 (a seeded draw from the Metagraph 100-studies table with the LOCALE training and validation runs excluded) plus the five genotype-B runs. All exist in Logan.
2. **Queries.** The first 1000 spots of each of the ten genotype-D runs were streamed from SRA with `fastq-dump` (mate 1 for paired runs); 500 were kept by a seeded shuffle. Stored sequences are the raw reads exactly as sequenced.
3. **Mutated copies.** The same queries with sequencing noise injected by [mutation-simulator](https://github.com/mkpython3/Mutation-Simulator) at SNP rate r with insertion and deletion rates r/10 each, r = 0.05 and 0.10. Every row keeps its `query_id` and ground-truth columns. The simulator's fasta and VCF per rate are under `mutations/`.

## Files

| File | Contents |
|---|---|
| `accs.txt` | the 55 indexed SRA run accessions, one per line (50 distractors, then the 5 relevant runs) |
| `queries.parquet` | one row per query with its ground truth (schema below) |
| `queries_mut0.00.parquet` | `queries.parquet` plus a `mutation_rate` column (clean) |
| `queries_mut0.05.parquet`, `queries_mut0.10.parquet` | same rows in the same order with `query_sequence` mutated |
| `queries.fa`, `queries.fa.fai` | the clean queries as fasta and its index |
| `mutations/` | mutation-simulator fasta and VCF for each rate |

## Schema of `queries.parquet`

| Column | Type | Meaning |
|---|---|---|
| `accession` | str | genotype-D run the read was sampled from (not indexed) |
| `query_id` | str | `<accession>.<spot>` from fastq-dump |
| `query_sequence` | str | the raw read |
| `contig_accession` | list[str] | **relevant set**: the five genotype-B runs, for every query |
| `contig_id` | list[str] | the relevant run accession (no contig-level ground truth) |
| `identity`, `ratio`, `contig_len`, `aln_interval_contig` | list | null |
| `strand` | str | null |

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
