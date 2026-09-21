"""Build a cross-genotype viral benchmark bundle with metadata-defined relevance.

The index set is an accession list of distractors (e.g. accessions/sra50.txt)
plus a small list of relevant runs of one genotype; queries are raw reads
sampled from runs of a different genotype. Every query's relevant set is the
full relevant list (Rq = its size): relevance comes from SRA metadata, no
aligner is involved. The bundle has the standard layout and the same schema
as the alignment bundles, so run_benchmark.py runs on it unchanged: identity,
ratio, contig_len, aln_interval_contig and strand are null, and contig_id
holds the relevant accession itself.

Usage:
    uv run python benchmark/build_viral_bundle.py \
        --index_accs accessions/sra50.txt \
        --relevant_accs accessions/hbv_genotype_B.txt \
        --query_accs accessions/hbv_genotype_D.txt \
        --work_dir constructed/sra55viral/work --output_path constructed/sra55viral/bundle \
        --seed=55
"""

from pathlib import Path

import polars as pl
from Bio import SeqIO
from finalize_query_dataset import write_output_fasta
from jsonargparse import auto_cli
from jsonargparse.typing import Path_fr
from mutate_queries import DEFAULT_RATES
from mutate_queries import main as mutate_bundle
from sample_raw_read_queries import main as sample_raw_reads

INTERVAL_DTYPE = pl.Struct(
    {"aln_start_index_contig": pl.Int64, "aln_end_index_contig": pl.Int64}
)


def read_list(path) -> list[str]:
    return Path(path).read_text().split()


def main(
    index_accs: Path_fr,
    relevant_accs: Path_fr,
    query_accs: Path_fr,
    work_dir: str,
    output_path: str,
    num_reads: int = 1000,
    queries_per_accession: int | None = None,
    final_query_num: int = 500,
    seed: int = 0,
    rates: list[float] = DEFAULT_RATES,
    num_download_procs: int = 4,
):
    """
    Args:
        index_accs: distractor accessions (e.g. accessions/sra50.txt).
        relevant_accs: accessions relevant to every query; indexed too.
        query_accs: runs to sample query reads from; not indexed.
        work_dir: raw reads and the initial query fasta are written here.
        output_path: bundle directory (created).
    """
    work_dir = Path(work_dir).resolve()
    output_path = Path(output_path).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    index_list = read_list(index_accs)
    relevant_list = read_list(relevant_accs)
    query_list = read_list(query_accs)
    all_accs = index_list + relevant_list
    assert len(all_accs) == len(set(all_accs)), "index/relevant lists overlap"
    assert not set(query_list) & set(all_accs), "query runs must not be indexed"

    queries_fa = work_dir / "queries.fa"
    sample_raw_reads(
        str(query_accs),
        str(work_dir / "raw_reads"),
        str(queries_fa),
        num_reads=num_reads,
        queries_per_accession=queries_per_accession,
        num_download_procs=num_download_procs,
        seed=seed,
    )

    queries_df = pl.from_dicts(
        [
            {
                "accession": record.id.split(".")[0],
                "query_id": record.id,
                "query_sequence": str(record.seq),
            }
            for record in SeqIO.parse(queries_fa, "fasta")
        ],
        schema={"accession": pl.Utf8, "query_id": pl.Utf8, "query_sequence": pl.Utf8},
    )
    assert final_query_num <= len(queries_df), (
        f"Requested {final_query_num} queries, only {len(queries_df)} sampled"
    )
    queries_df = queries_df.sort("query_id").sample(
        final_query_num, shuffle=True, seed=seed
    )

    k = len(relevant_list)
    nulls = [None] * k
    queries_df = queries_df.with_columns(
        pl.lit(relevant_list, dtype=pl.List(pl.Utf8)).alias("contig_id"),
        pl.lit(relevant_list, dtype=pl.List(pl.Utf8)).alias("contig_accession"),
        pl.lit(nulls, dtype=pl.List(pl.Float64)).alias("identity"),
        pl.lit(nulls, dtype=pl.List(pl.Float64)).alias("ratio"),
        pl.lit(nulls, dtype=pl.List(pl.Int64)).alias("contig_len"),
        pl.lit(nulls, dtype=pl.List(INTERVAL_DTYPE)).alias("aln_interval_contig"),
        pl.lit(None, dtype=pl.Utf8).alias("strand"),
    )

    output_path.mkdir(parents=True, exist_ok=True)
    (output_path / "accs.txt").write_text("".join(a + "\n" for a in all_accs))
    queries_df.write_parquet(output_path / "queries.parquet")
    write_output_fasta(queries_df, output_path / "queries.fa")
    print(f"Wrote {len(queries_df)} queries with Rq={k} to {output_path}")

    mutate_bundle(str(output_path), rates=rates)


if __name__ == "__main__":
    auto_cli(main, as_positional=False)
