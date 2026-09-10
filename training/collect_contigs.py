import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import polars as pl
from Bio import SeqIO
from jsonargparse import auto_cli
from tqdm import tqdm


def collect_and_filter_contigs(
    contig_file: Path, min_contig_len: int, max_contig_len: int
):
    srr_id = contig_file.parent.stem
    seqs = []
    for record in SeqIO.parse(contig_file, "fasta"):
        seq = str(record.seq)
        if len(seq) >= max_contig_len:
            chunks = [
                seq[i : i + max_contig_len] for i in range(0, len(seq), max_contig_len)
            ]
            seqs.extend(c for c in chunks if len(c) >= min_contig_len)
        elif len(seq) >= min_contig_len:
            seqs.append(seq)

    return seqs, srr_id


def main(
    contig_dir: str,
    output_path: str,
    min_contig_len: int = 200,
    max_contig_len: int = 1024,
):
    """Collect downloaded logan contigs into a training parquet.

    Reads every <acc>/<acc>.contigs.fa under contig_dir (the output of
    download_logan_contigs.py), splits contigs longer than max_contig_len into
    max_contig_len chunks, drops anything shorter than min_contig_len, shuffles,
    and writes a parquet with accession/sequence/sequence_len columns — the
    format the locale Batcher consumes in "contig" mode.
    """
    os.environ["POLARS_MAX_THREADS"] = "1"
    contig_dir: Path = Path(contig_dir).resolve()
    output_path: Path = Path(output_path).resolve()
    contig_files = list(contig_dir.rglob("*.contigs.fa"))
    all_contigs = []
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = [
            executor.submit(
                collect_and_filter_contigs, contig_file, min_contig_len, max_contig_len
            )
            for contig_file in contig_files
        ]
        for future in tqdm(
            as_completed(futures),
            total=len(contig_files),
            desc="Collecting contig files",
        ):
            seqs, srr_id = future.result()
            all_contigs.append(pl.DataFrame({"accession": srr_id, "sequence": seqs}))

    all_contigs = pl.concat(all_contigs)
    print("Shuffling...")
    all_contigs = all_contigs.sample(fraction=1.0, shuffle=True, seed=42)
    all_contigs = all_contigs.with_columns(
        pl.col("sequence").str.len_chars().alias("sequence_len")
    )
    print(f"Writing {len(all_contigs)} contigs to {output_path}")
    all_contigs.write_parquet(output_path)


if __name__ == "__main__":
    auto_cli(main)
