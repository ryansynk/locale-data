"""Summarise a reference-genome data package and/or a cropped parquet.

With an assembly report: species / accession / total bp per genome as LaTeX
(the paper's reference-genome table). With a parquet: crops per organism.
"""

import json
from pathlib import Path

import polars as pl
from jsonargparse import auto_cli


def main(
    assembly_data_report_path: Path | None = None,
    parquet_path: Path | None = None,
):
    if assembly_data_report_path is not None:
        with open(assembly_data_report_path) as f:
            records = [json.loads(line) for line in f]
        rows = [
            {
                "species": r["organism"]["organismName"],
                "accession": r["accession"],
                "total_bp": r["assemblyStats"]["totalSequenceLength"],
            }
            for r in records
        ]
        print("\\begin{tabular}{llr}\n\\toprule\nspecies & accession & total\\_bp \\\\\n\\midrule")
        for r in rows:
            print(f"{r['species']} & {r['accession']} & {r['total_bp']} \\\\")
        print("\\bottomrule\n\\end{tabular}")

    if parquet_path is not None:
        df = pl.scan_parquet(parquet_path)
        stats = (
            df.group_by("accession", "organism")
            .agg(pl.len().alias("crops"), pl.col("sequence").str.len_chars().sum().alias("bp"))
            .sort("crops")
            .collect()
        )
        print(stats)
        print(f"total crops: {stats['crops'].sum():,}   total bp: {stats['bp'].sum():,}")


if __name__ == "__main__":
    auto_cli(main)
