"""Write one pre-mutated query file per mutation rate into a benchmark bundle.

Runs after finalize_query_dataset.py. For every rate r, mutation-simulator is
run on queries.fa with SNP rate r and insertion and deletion rates r/10 each,
and the result is written as queries_mut<r>.parquet: queries.parquet with
query_sequence replaced by the mutated read and a mutation_rate column added.
Row order and query_id are identical to queries.parquet, so the benchmark's
seeded subsample selects the same queries at every rate. Rate 0 is written
without running the simulator.

mutation-simulator has no seed flag, so the mutated fasta and VCF it emits are
kept under <bundle>/mutations/ as the provenance record; regenerating a rate
produces different mutations, and the parquet is what the benchmark reads.

Usage:
    uv run python mutate_queries.py --bundle <bundle_dir>
"""

import subprocess
from pathlib import Path

import polars as pl
from Bio import SeqIO
from jsonargparse import auto_cli
from jsonargparse.typing import Path_drw

DEFAULT_RATES = [0.0, 0.05, 0.10]


def rate_tag(rate: float) -> str:
    return f"mut{rate:.2f}"


def run_mutation_simulator(
    queries_fa: Path, out_dir: Path, rate: float, indel_fraction: float
) -> Path:
    """Mutate queries_fa at the given SNP rate; return the mutated fasta path.

    mutation-simulator truncates the output basename at its first '.', so the
    basename is dot-free (mut05pct, not mut0.05). Indel lengths use the
    simulator defaults (1 to 2 bases).
    """
    out_base = out_dir / f"mut{round(rate * 100):02d}pct"
    cmd = [
        "mutation-simulator",
        "--no-progress",
        "--no-color",
        "--quiet",
        "--output",
        str(out_base),
        str(queries_fa),
        "args",
        "--snp",
        str(rate),
        "--insert",
        str(rate * indel_fraction),
        "--deletion",
        str(rate * indel_fraction),
    ]
    subprocess.run(cmd, check=True)
    mutated_fa = out_base.with_name(out_base.name + "_ms.fa")
    if not mutated_fa.exists():
        raise FileNotFoundError(f"mutation-simulator wrote no fasta at {mutated_fa}")
    return mutated_fa


def read_sequences(fasta: Path) -> pl.DataFrame:
    ids, seqs = [], []
    for record in SeqIO.parse(fasta, "fasta"):
        ids.append(record.id)
        seqs.append(str(record.seq))
    return pl.DataFrame({"query_id": ids, "query_sequence": seqs})


def mutated_queries(
    queries: pl.DataFrame, sequences: pl.DataFrame, rate: float
) -> pl.DataFrame:
    """queries with query_sequence swapped for sequences', in queries' row order."""
    out = queries.drop("query_sequence").join(
        sequences, on="query_id", how="left", maintain_order="left"
    )
    if out["query_id"].to_list() != queries["query_id"].to_list():
        raise ValueError("query_id order changed during join")
    if out["query_sequence"].null_count():
        raise ValueError("some queries have no mutated sequence")
    return out.select(queries.columns).with_columns(
        pl.lit(rate, dtype=pl.Float64).alias("mutation_rate")
    )


def main(
    bundle: Path_drw,
    rates: list[float] = DEFAULT_RATES,
    indel_fraction: float = 0.1,
):
    """
    Args:
        bundle: dataset directory holding queries.parquet and queries.fa.
        rates: SNP rates; one queries_mut<rate>.parquet is written per rate.
        indel_fraction: insertion and deletion rates are each this fraction
            of the SNP rate.
    """
    bundle = Path(bundle)
    queries = pl.read_parquet(bundle / "queries.parquet")
    queries_fa = bundle / "queries.fa"
    mutations_dir = bundle / "mutations"
    mutations_dir.mkdir(exist_ok=True)

    for rate in rates:
        if rate == 0.0:
            sequences = queries.select("query_id", "query_sequence")
        else:
            sequences = read_sequences(
                run_mutation_simulator(queries_fa, mutations_dir, rate, indel_fraction)
            )
        out = mutated_queries(queries, sequences, rate)
        out_path = bundle / f"queries_{rate_tag(rate)}.parquet"
        out.write_parquet(out_path)
        changed = (out["query_sequence"] != queries["query_sequence"]).sum()
        print(f"{rate_tag(rate)}: {len(out)} queries, {changed} changed -> {out_path}")


if __name__ == "__main__":
    auto_cli(main, as_positional=False)
