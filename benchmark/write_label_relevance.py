"""Write a label-defined relevance file for a bundle's queries.

Every query gets the same relevant set: the accessions in `relevant_accs`
(e.g. all HBV runs of a cross-genotype bundle), regardless of alignment. The
output has the columns the benchmark's scorer reads from queries.parquet
(query_id, contig_accession, identity), so it can be passed to
scripts/print_results.py in place of the bundle's alignment-defined queries
file to score the same search results under label relevance:

    uv run python benchmark/write_label_relevance.py \
        --bundle constructed/sra64viral-v2/bundle \
        --relevant_accs "[accessions/hbv_genotype_A.txt,accessions/hbv_genotype_B.txt,accessions/hbv_genotype_C.txt]"
    # then, in locale/benchmark:
    uv run python scripts/print_results.py results/sra64viral \
        <bundle>/label_relevance.parquet <bundle>/accs.txt
"""

from pathlib import Path

import polars as pl
from jsonargparse import auto_cli


def main(bundle: str, relevant_accs: list[str], output: str | None = None):
    """
    Args:
        bundle: bundle directory holding queries.parquet and accs.txt.
        relevant_accs: accession list files; their union is relevant to
            every query. Each accession must be indexed (in accs.txt).
        output: defaults to <bundle>/label_relevance.parquet.
    """
    bundle = Path(bundle)
    relevant = [a for f in relevant_accs for a in Path(f).read_text().split()]
    assert len(relevant) == len(set(relevant)), "duplicate relevant accessions"
    indexed = set((bundle / "accs.txt").read_text().split())
    assert set(relevant) <= indexed, sorted(set(relevant) - indexed)
    queries = pl.read_parquet(bundle / "queries.parquet", columns=["query_id"])
    out = queries.with_columns(
        pl.lit(relevant, dtype=pl.List(pl.Utf8)).alias("contig_accession"),
        pl.lit([1.0] * len(relevant), dtype=pl.List(pl.Float64)).alias("identity"),
    )
    output = Path(output) if output else bundle / "label_relevance.parquet"
    out.write_parquet(output)
    print(
        f"Wrote {len(out)} queries with {len(relevant)} relevant accessions to {output}"
    )


if __name__ == "__main__":
    auto_cli(main, as_positional=False)
