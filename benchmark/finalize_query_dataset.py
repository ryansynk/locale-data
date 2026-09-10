from pathlib import Path
import math
import shutil
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


def has_overlap(intervals, chunk_size, mode):
    """Helper function for filtering out queries that align to similar regions on the same contig.

    Given a list of alignment intervals on a contig, determines whether any pair
    of intervals are close enough to be considered overlapping, based on the
    filtering mode.

    Args:
        intervals: List of structs with two int fields representing
            alignment start and end positions on a contig.
        chunk_size: The size of the chunk window used for the strict
            overlap check. Default is 1024.
        mode: Filtering stringency for overlap detection.
            - "strict": Flags overlap if the gap between any two intervals
              is less than chunk_size, assuming worst-case chunk placement.
            - "lenient": Flags overlap only when intervals directly overlap.
            - "none": No overlap filtering; always returns False.

    Returns:
        bool: True if any intervals are considered overlapping under the
        given mode, False otherwise.
    """
    if len(intervals) <= 1:
        return False

    sorted_ivs = sorted(
        (iv["aln_start_index_contig"], iv["aln_end_index_contig"]) for iv in intervals
    )

    for i in range(1, len(sorted_ivs)):
        gap = sorted_ivs[i][0] - sorted_ivs[i - 1][1]

        if mode == "strict":
            # Worst case chunk placement — any gap < chunk_size is suspect
            if gap < chunk_size:
                return True
        elif mode == "lenient":
            # Only flag direct overlap
            if gap < 0:
                return True
        # mode == "none" — no overlap check, always False

    return False


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
    stats_df = pl.read_parquet(minimap_stats)
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
    stats_df = stats_df.join(
        queries_df.select(
            [
                pl.col("read_id"),
                pl.col("query_sequence").str.len_chars().alias("read_length"),
            ]
        ),
        on="read_id",
    )
    stats_df = stats_df.with_columns(
        (pl.col("matches") / pl.col("read_length")).alias("identity")
    )
    stats_df = stats_df.with_columns(
        (pl.col("target_len") / pl.col("read_length")).alias("ratio")
    )
    stats_df = stats_df.filter(pl.col("identity") > identity_filter)
    stats_df = stats_df.with_columns(
        pl.struct("aln_start_index_contig", "aln_end_index_contig").alias(
            "aln_interval_contig"
        )
    )
    stats_df = stats_df.group_by("read_id").agg(
        [
            "contig_id",
            "accession",
            "identity",
            "ratio",
            "target_len",
            "aln_interval_contig",
        ]
    )
    max_matches = math.ceil(num_indexed_accessions * max_match_percent)
    stats_df = stats_df.filter(pl.col("accession").list.len() <= max_matches)
    stats_df = stats_df.rename(
        {"accession": "contig_accession", "target_len": "contig_len"}
    )
    queries_df = queries_df.join(stats_df, on="read_id")
    queries_df = queries_df.rename({"read_id": "query_id"})

    if filter_mode != "none":
        overlapping_queries = (
            queries_df.explode("contig_id", "aln_interval_contig")
            .group_by("contig_id")
            .agg(pl.col("query_id"), pl.col("aln_interval_contig"))
            .with_columns(
                pl.col("aln_interval_contig")
                .map_elements(
                    lambda intervals: has_overlap(
                        intervals, chunk_size=chunk_size, mode=filter_mode
                    ),
                    return_dtype=pl.Boolean,
                )
                .alias("has_overlap")
            )
            .filter(pl.col("has_overlap"))
            .explode("query_id")["query_id"]
            .unique()
            .to_list()
        )
        queries_df = queries_df.filter(~pl.col("query_id").is_in(overlapping_queries))

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


if __name__ == "__main__":
    auto_cli(main)
