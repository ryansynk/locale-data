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
    seed: int = 0,
):
    """Filter query alignments into a final benchmark dataset bundle.

    Writes a self-contained dataset directory at output_path containing
    queries.parquet, queries.fa, and a copy of the indexed accession list as
    accs.txt — the exact layout the benchmark consumes.

    Alignments with identity (exact matches / read length) <= identity_filter
    are dropped first. Then, per query, in this order:
      1. drop it if no surviving alignment is to the accession it was sampled
         from (its source must be in its relevant set);
      2. drop it if it aligns to more than ceil(max_match_percent * number of
         indexed accessions) accessions;
      3. drop it if it overlaps another query (filter_mode, see
         overlapping_query_ids).
    The survivors are shuffled with `seed` and the first final_query_num kept,
    so the same inputs and seed reproduce queries.parquet exactly.

    The stored query is the raw read; alignment ran on both strands and the
    `strand` column is the strand of the best hit to the source accession.
    """
    minimap_stats: Path = Path(minimap_stats)
    queries_fa: Path = Path(queries_fa)
    accessions_file: Path = Path(accessions_file)
    accessions = accessions_file.read_text().split()
    assert len(accessions) == len(set(accessions)), "duplicate accessions in list"
    num_indexed_accessions = len(accessions)
    if "strand" not in pl.scan_parquet(minimap_stats).collect_schema():
        raise SystemExit(
            f"{minimap_stats} has no `strand` column: it was produced with the "
            "forward-only aligner. Re-run get_query_alignments_minimap.py."
        )

    queries_rows = []
    for record in SeqIO.parse(queries_fa, "fasta"):
        queries_rows.append(
            {
                "accession": str(record.id.split(".")[0]),
                "read_id": str(record.id),
                "query_sequence": str(record.seq),
            }
        )
    queries_df = pl.from_dicts(
        queries_rows,
        schema={"accession": pl.Utf8, "read_id": pl.Utf8, "query_sequence": pl.Utf8},
    )
    initial = len(queries_df)

    read_info = queries_df.select(
        pl.col("read_id"),
        pl.col("accession").alias("source_accession"),
        pl.col("query_sequence").str.len_chars().alias("read_length"),
    ).lazy()

    # Lazy pipeline over the (very large) stats parquet: only needed columns
    # are read, and filters run before anything is fully materialized. Rows
    # are sorted so the per-query lists come out in a fixed order.
    is_source = pl.col("accession") == pl.col("source_accession")
    stats_df = (
        pl.scan_parquet(minimap_stats)
        .join(read_info, on="read_id")
        .with_columns(
            (pl.col("matches") / pl.col("read_length")).alias("identity"),
            (pl.col("target_len") / pl.col("read_length")).alias("ratio"),
        )
        .filter(pl.col("identity") > identity_filter)
        .sort("read_id", "accession", "contig_id", "aln_start_index_contig")
        .with_columns(
            pl.struct("aln_start_index_contig", "aln_end_index_contig").alias(
                "aln_interval_contig"
            )
        )
        .group_by("read_id", maintain_order=True)
        .agg(
            [
                "contig_id",
                "accession",
                "identity",
                "ratio",
                "target_len",
                "aln_interval_contig",
                pl.col("strand").filter(is_source).first().alias("strand"),
            ]
        )
        .rename({"accession": "contig_accession", "target_len": "contig_len"})
        .collect()
    )
    queries_df = queries_df.join(stats_df, on="read_id")
    queries_df = queries_df.rename({"read_id": "query_id"})

    # 1. The source accession must be in the relevant set
    queries_df = queries_df.filter(
        pl.col("contig_accession").list.contains(pl.col("accession"))
    )
    n_with_source = len(queries_df)
    # 2. Max-match
    max_matches = math.ceil(num_indexed_accessions * max_match_percent)
    queries_df = queries_df.filter(
        pl.col("contig_accession").list.len() <= max_matches
    )
    n_within_max = len(queries_df)
    # 3. Overlap
    if filter_mode != "none":
        overlapping = overlapping_query_ids(
            queries_df, chunk_size=chunk_size, mode=filter_mode
        )
        queries_df = queries_df.join(overlapping, on="query_id", how="anti")
    eligible = len(queries_df)

    assert final_query_num <= eligible, (
        f"Requested {final_query_num} queries, only {eligible} after filtering"
    )
    queries_df = queries_df.sort("query_id").sample(
        final_query_num, shuffle=True, seed=seed
    )
    assert queries_df["strand"].null_count() == 0

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
    n_minus = (queries_df["strand"] == "-").sum()
    print(
        f"Wrote {len(queries_df)} queries to {output_path}: {initial} initial, "
        f"-{initial - n_with_source} no source hit, "
        f"-{n_with_source - n_within_max} max-match (>{max_matches}), "
        f"-{n_within_max - eligible} overlap, {eligible} eligible; "
        f"{n_minus}/{len(queries_df)} on the - strand"
    )


if __name__ == "__main__":
    auto_cli(main)
