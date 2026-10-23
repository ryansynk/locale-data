import random
import shutil
from pathlib import Path

import polars as pl
from Bio import SeqIO
from jsonargparse import auto_cli


def main(
    input_dataset: Path,
    output_dataset: Path,
    hit_contigs_path: Path,
    subsample_fraction: float = 0.1,
    seed: int = 0,
):
    # First copy over all files other than the logan_accessions_dir
    for f in input_dataset.iterdir():
        if f.is_file():
            shutil.copy(f, output_dataset / f.name)

    # Collect all contigs needed to be sampled (positives and hard negatives)
    truth_contigs = (
        pl.read_parquet(input_dataset / "queries.parquet")
        .select(pl.col("contig_accession").alias("accession"), "contig_id")
        .explode("accession", "contig_id")
        .unique()
    )  # (accession, contig_id), same shape as hit_contigs
    hit_contigs = pl.read_parquet(hit_contigs_path)  # (accession, contig_id)
    matches = (
        pl.concat([truth_contigs, hit_contigs]).group_by("accession").agg("contig_id")
    )

    # stream over all fastas, collected needed contigs and subsample otherwise
    keep = {acc: set(contig_ids) for acc, contig_ids in matches.iter_rows()}
    rng = random.Random(seed)
    for fasta in sorted((input_dataset / "logan_accessions").glob("*/*.fa")):
        acc = fasta.parent.name
        wanted = keep.get(acc, set())
        out = output_dataset / "logan_accessions" / acc / fasta.name
        out.parent.mkdir(parents=True, exist_ok=True)
        n_in = n_out = 0
        with open(out, "w") as f:
            for rec in SeqIO.parse(fasta, "fasta"):
                n_in += 1
                if rec.id in wanted or rng.random() < subsample_fraction:
                    SeqIO.write(rec, f, "fasta")
                    n_out += 1
        print(f"{acc}: kept {n_out}/{n_in} contigs ({len(wanted)} wanted)")


if __name__ == "__main__":
    auto_cli(main, as_positional=False)
