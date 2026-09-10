import os
from typing import Literal
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import polars as pl
import pysam
from jsonargparse import auto_cli
from jsonargparse.typing import Path_drw, Path_fr, Path_dc
from tqdm import tqdm


def run_minimap(queries_fa: str, contig_fa: str, output_sam: str):
    cmd = ["minimap2", "-a", "--eqx", "--for-only", contig_fa, queries_fa]
    with open(output_sam, "w") as out_file:
        subprocess.run(cmd, stdout=out_file, stderr=subprocess.DEVNULL, check=True)


def get_best_alignments(input_sam):
    accession = input_sam.name.split(".")[0]
    best_alns = {}

    # Open the SAM file for reading
    with pysam.AlignmentFile(input_sam, "r") as infile:
        # fetch(until_eof=True) allows iterating over SAM files without an index
        for read in infile.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            # Reconstruct the full read length from the CIGAR string
            query_len = read.infer_query_length()
            if not query_len:
                continue

            # get_cigar_stats()[0] returns an array of total operation counts:
            # Index 0: M, Index 1: I, Index 2: D, Index 7: = (exact match)
            cigar_stats = read.get_cigar_stats()[0]
            eq_bases = cigar_stats[7]
            assert eq_bases > 0
            exact_matches = eq_bases

            # Calculate Identity
            ref_name = read.reference_name
            ref_length = infile.lengths[read.reference_id]
            is_primary = not (read.is_secondary or read.is_supplementary)
            qname = read.query_name

            # Keep track of the highest identity read per query
            current_best = best_alns.get(qname)

            aln_start = read.reference_start
            aln_end = read.reference_end

            if current_best is None:
                best_alns[qname] = {
                    "read_id": qname,
                    "is_primary": is_primary,
                    "matches": exact_matches,
                    "accession": accession,
                    "contig_id": ref_name,
                    "target_len": ref_length,
                    "aln_start_index_contig": aln_start,
                    "aln_end_index_contig": aln_end,
                }
            elif exact_matches > current_best["matches"]:
                best_alns[qname] = {
                    "read_id": qname,
                    "is_primary": is_primary,
                    "matches": exact_matches,
                    "accession": accession,
                    "contig_id": ref_name,
                    "target_len": ref_length,
                    "aln_start_index_contig": aln_start,
                    "aln_end_index_contig": aln_end,
                }
            elif exact_matches == current_best["matches"]:
                # Tie-breaker: If identities are identical, prefer the primary alignment flag
                if is_primary and not current_best["is_primary"]:
                    best_alns[qname] = {
                        "read_id": qname,
                        "is_primary": is_primary,
                        "matches": exact_matches,
                        "accession": accession,
                        "contig_id": ref_name,
                        "target_len": ref_length,
                        "aln_start_index_contig": aln_start,
                        "aln_end_index_contig": aln_end,
                    }

    return list(best_alns.values())


def main(
    queries_fa: Path_fr,
    logan_path: Path_drw,
    output_path: Path_dc,
    num_workers: int,
    chunk_size: int = 10_000,
):
    os.environ["POLARS_MAX_THREADS"] = "1"
    logan_path: Path = Path(logan_path).resolve()
    output_path: Path = Path(output_path).resolve()
    queries_fa: Path = Path(queries_fa).resolve()
    # queries_df = pl.read_parquet(queries_parquet)
    contig_files = list(logan_path.rglob("*.contigs.fa"))

    alignments_dir = Path(output_path / "alignments")
    alignments_dir.mkdir(exist_ok=True)
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for contig_file in contig_files:
            accession = contig_file.parent.stem
            sam_path = alignments_dir / f"{accession}.sam"
            futures.append(
                executor.submit(
                    run_minimap, str(queries_fa), str(contig_file), str(sam_path)
                )
            )

        for f in tqdm(
            as_completed(futures), total=len(futures), desc="Running minimap2..."
        ):
            f.result()

    aln_files = list(alignments_dir.glob("*.sam"))
    best_alignments = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        futures = [
            executor.submit(get_best_alignments, aln_file) for aln_file in aln_files
        ]
        for f in tqdm(
            as_completed(futures), total=len(futures), desc="Extracting best alignments"
        ):
            best_alignments.extend(f.result())
    # best_alignments = []
    # for aln_file in tqdm(aln_files, desc="Extracting best alignments..."):
    #    best_alignments.extend(get_best_alignments(aln_file))

    parquet_out = (
        output_path / f"{queries_fa.name.split('.')[0]}_alignment_stats.parquet"
    )
    best_alignments = pl.from_dicts(best_alignments)
    best_alignments.write_parquet(parquet_out)


if __name__ == "__main__":
    auto_cli(main)
