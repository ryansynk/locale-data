import math
import shutil
from pathlib import Path
from typing import Literal

import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from jsonargparse import auto_cli
from jsonargparse.typing import Path_dc, Path_fr


def write_output_fasta(queries_df: pl.DataFrame, output_fasta: Path):
    records = (
        SeqRecord(Seq(sequence), id=read_id, description="")
        for read_id, sequence in queries_df.select(
            ["query_id", "query_sequence"]
        ).iter_rows()
    )
    SeqIO.write(records, output_fasta, "fasta")


def overlapping_query_ids(
    queries_df: pl.DataFrame, chunk_size: int, mode: str
) -> pl.DataFrame:
    """Find queries that align to similar regions on the same contig.

    A contig is flagged when any two of its alignment intervals (sorted by
    start) are close enough to be considered overlapping:
      - "strict": the gap between consecutive intervals is < chunk_size,
        assuming worst-case chunk placement.
      - "lenient": intervals directly overlap (negative gap).
    Every query with an alignment on a flagged contig is returned.

    Only the columns needed for the check are exploded: exploding the full
    frame would duplicate the query sequence and every other list column per
    alignment row, which does not fit in memory at this scale.
    """
    intervals = (
        queries_df.lazy()
        .select("query_id", "contig_id", "aln_interval_contig")
        .explode("contig_id", "aln_interval_contig")
        .unnest("aln_interval_contig")
        .sort("contig_id", "aln_start_index_contig", "aln_end_index_contig")
    )
    gap_threshold = chunk_size if mode == "strict" else 0
    flagged_contigs = (
        intervals.with_columns(
            (
                pl.col("aln_start_index_contig")
                - pl.col("aln_end_index_contig").shift(1).over("contig_id")
            ).alias("gap")
        )
        .filter(pl.col("gap") < gap_threshold)
        .select("contig_id")
        .unique()
    )
    return (
        intervals.join(flagged_contigs, on="contig_id")
        .select("query_id")
        .unique()
        .collect()
    )


def main(
    queries_fa: Path_fr,
    minimap_stats: Path_fr,
    accessions_file: Path_fr,
    output_path: Path_dc,
    chunk_size: int = 1024,
    filter_mode: Literal["strict", "lenient", "none"] = "strict",
    identity_filter: float = 0.9,
    max_match_percent: float = 0.2,
    final_query_num: int = 500,
):
    """Filter query alignments into a final benchmark dataset bundle.

    Writes a self-contained dataset directory at output_path containing
    queries.parquet, queries.fa, and a copy of the indexed accession list as
    accs.txt — the exact layout the benchmark consumes.
    """
    minimap_stats: Path = Path(minimap_stats)
    queries_fa: Path = Path(queries_fa)
    accessions_file: Path = Path(accessions_file)
    accessions = accessions_file.read_text().splitlines()
    num_indexed_accessions = len(accessions)

    queries_rows = []
    for record in SeqIO.parse(queries_fa, "fasta"):
        queries_rows.append(
            {
                "accession": str(record.id.split(".")[0]),
                "read_id": str(record.id),
                "query_sequence": str(record.seq),
            }
        )
    queries_df = pl.from_dicts(queries_rows)

    read_lengths = queries_df.select(
        pl.col("read_id"),
        pl.col("query_sequence").str.len_chars().alias("read_length"),
    ).lazy()

    # Lazy pipeline over the (very large) stats parquet: only needed columns
    # are read, and filters run before anything is fully materialized
    max_matches = math.ceil(num_indexed_accessions * max_match_percent)
    stats_df = (
        pl.scan_parquet(minimap_stats)
        .join(read_lengths, on="read_id")
        .with_columns(
            (pl.col("matches") / pl.col("read_length")).alias("identity"),
            (pl.col("target_len") / pl.col("read_length")).alias("ratio"),
        )
        .filter(pl.col("identity") > identity_filter)
        .with_columns(
            pl.struct("aln_start_index_contig", "aln_end_index_contig").alias(
                "aln_interval_contig"
            )
        )
        .group_by("read_id")
        .agg(
            [
                "contig_id",
                "accession",
                "identity",
                "ratio",
                "target_len",
                "aln_interval_contig",
            ]
        )
        .filter(pl.col("accession").list.len() <= max_matches)
        .rename({"accession": "contig_accession", "target_len": "contig_len"})
        .collect()
    )
    queries_df = queries_df.join(stats_df, on="read_id")
    queries_df = queries_df.rename({"read_id": "query_id"})

    if filter_mode != "none":
        overlapping = overlapping_query_ids(
            queries_df, chunk_size=chunk_size, mode=filter_mode
        )
        queries_df = queries_df.join(overlapping, on="query_id", how="anti")

    # Filter queries not aligned to their original accession
    queries_df = queries_df.filter(
        ~pl.col("contig_accession").list.contains(pl.col("accession"))
    )
    assert final_query_num <= len(queries_df), (
        f"Requested {final_query_num} queries, only {len(queries_df)} after filtering"
    )
    queries_df = queries_df.sample(final_query_num, shuffle=True)

    # Every accession the final queries reference must come from the indexed set
    acc_set = set(accessions)
    referenced_accs = set(queries_df["accession"].to_list()) | set(
        queries_df["contig_accession"].explode().to_list()
    )
    assert referenced_accs <= acc_set, (
        f"Queries reference accessions not in {accessions_file}: "
        f"{sorted(referenced_accs - acc_set)[:5]}"
    )

    output_path: Path = Path(output_path).resolve()
    output_path.mkdir(parents=True, exist_ok=True)
    queries_df.write_parquet(output_path / "queries.parquet")
    write_output_fasta(queries_df, output_path / "queries.fa")
    shutil.copy(accessions_file, output_path / "accs.txt")
    print(f"Wrote {len(queries_df)} queries to {output_path}")


if __name__ == "__main__":
    auto_cli(main)
