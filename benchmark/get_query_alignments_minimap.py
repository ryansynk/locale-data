import os
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Must be set before polars is imported to actually cap its thread pool;
# worker processes fork from this state, so this keeps them single-threaded.
os.environ["POLARS_MAX_THREADS"] = "1"

import polars as pl
from jsonargparse import auto_cli
from jsonargparse.typing import Path_dc, Path_drw, Path_fr
from tqdm import tqdm

CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")

BEST_ALN_SCHEMA = {
    "read_id": pl.Utf8,
    "is_primary": pl.Boolean,
    "matches": pl.Int64,
    "accession": pl.Utf8,
    "contig_id": pl.Utf8,
    "target_len": pl.Int64,
    "aln_start_index_contig": pl.Int64,
    "aln_end_index_contig": pl.Int64,
}


def run_minimap(queries_fa: str, contig_fa: str, output_sam: str):
    # A SAM that already exists is complete (partial output goes to .tmp),
    # so reruns reuse prior work instead of re-aligning.
    output_sam = Path(output_sam)
    if output_sam.exists():
        return
    # -t 1: parallelism comes from one process per contig file; -K 200M caps
    # the in-flight query batch so hundreds of workers fit in node memory.
    # --sam-hit-only: unmapped records are skipped downstream anyway and
    # would add ~1GB per SAM.
    cmd = [
        "minimap2",
        "-a",
        "--eqx",
        "--for-only",
        "--sam-hit-only",
        "-t",
        "1",
        "-K",
        "200M",
        contig_fa,
        queries_fa,
    ]
    tmp_sam = output_sam.with_suffix(".sam.tmp")
    with open(tmp_sam, "w") as out_file:
        subprocess.run(cmd, stdout=out_file, stderr=subprocess.DEVNULL, check=True)
    tmp_sam.rename(output_sam)


def get_best_alignments(input_sam: Path, output_parquet: Path):
    """Reduce one SAM to the best alignment per query and write it as parquet.

    Parses the SAM as text rather than with pysam: Logan assemblies can have
    >10M contigs, and pysam materializes every @SQ header line into Python
    objects on open (GBs of RAM and minutes of CPU per file). Only the @SQ
    lengths of contigs that actually won a best alignment are looked up, in a
    second cheap scan over the header.
    """
    accession = input_sam.name.split(".")[0]
    best_alns = {}

    with open(input_sam) as infile:
        for line in infile:
            if line.startswith("@"):
                continue
            fields = line.split("\t", 6)
            flag = int(fields[1])
            if flag & 0x4:  # unmapped
                continue

            # --eqx CIGARs: '=' ops are exact matches; M/D/N/=/X consume ref
            eq_bases = 0
            ref_span = 0
            for length, op in CIGAR_OP.findall(fields[5]):
                n = int(length)
                if op == "=":
                    eq_bases += n
                    ref_span += n
                elif op in "MDNX":
                    ref_span += n
            assert eq_bases > 0
            exact_matches = eq_bases

            qname = fields[0]
            is_primary = not (flag & 0x900)  # neither secondary nor supplementary
            aln_start = int(fields[3]) - 1  # SAM POS is 1-based
            aln_end = aln_start + ref_span

            # Keep the highest-exact-match alignment per query; on ties,
            # prefer the primary alignment
            current_best = best_alns.get(qname)
            if (
                current_best is None
                or exact_matches > current_best["matches"]
                or (
                    exact_matches == current_best["matches"]
                    and is_primary
                    and not current_best["is_primary"]
                )
            ):
                best_alns[qname] = {
                    "read_id": qname,
                    "is_primary": is_primary,
                    "matches": exact_matches,
                    "accession": accession,
                    "contig_id": fields[2],
                    "target_len": None,
                    "aln_start_index_contig": aln_start,
                    "aln_end_index_contig": aln_end,
                }

    # Second pass over just the header for the lengths of winning contigs
    needed = {v["contig_id"] for v in best_alns.values()}
    lengths = {}
    with open(input_sam) as infile:
        for line in infile:
            if not line.startswith("@"):
                break
            if not line.startswith("@SQ"):
                continue
            sn = ln = None
            for tag in line.rstrip("\n").split("\t")[1:]:
                if tag.startswith("SN:"):
                    sn = tag[3:]
                elif tag.startswith("LN:"):
                    ln = int(tag[3:])
            if sn in needed:
                lengths[sn] = ln
    for v in best_alns.values():
        v["target_len"] = lengths[v["contig_id"]]

    if best_alns:
        df = pl.from_dicts(list(best_alns.values()), schema=BEST_ALN_SCHEMA)
    else:
        df = pl.DataFrame(schema=BEST_ALN_SCHEMA)
    tmp_parquet = output_parquet.with_suffix(".parquet.tmp")
    df.write_parquet(tmp_parquet)
    tmp_parquet.rename(output_parquet)


def main(
    queries_fa: Path_fr,
    logan_path: Path_drw,
    output_path: Path_dc,
    num_workers: int,
    chunk_size: int = 10_000,
):
    logan_path: Path = Path(logan_path).resolve()
    output_path: Path = Path(output_path).resolve()
    queries_fa: Path = Path(queries_fa).resolve()
    # Largest files first so a huge contig file doesn't serialize the tail
    contig_files = sorted(
        logan_path.rglob("*.contigs.fa"), key=lambda p: p.stat().st_size, reverse=True
    )

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

    best_dir = output_path / "best_alignments"
    best_dir.mkdir(exist_ok=True)
    aln_files = sorted(
        alignments_dir.glob("*.sam"), key=lambda p: p.stat().st_size, reverse=True
    )
    # A per-accession parquet that already exists is complete, so reruns
    # only parse SAMs that haven't been reduced yet
    todo = [
        (sam, best_dir / f"{sam.stem}.parquet")
        for sam in aln_files
        if not (best_dir / f"{sam.stem}.parquet").exists()
    ]
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(get_best_alignments, sam, pq) for sam, pq in todo
        ]
        for f in tqdm(
            as_completed(futures), total=len(futures), desc="Extracting best alignments"
        ):
            f.result()

    parquet_out = (
        output_path / f"{queries_fa.name.split('.')[0]}_alignment_stats.parquet"
    )
    # Stream the per-accession parquets into one file without holding all
    # rows in memory at once
    pl.scan_parquet(best_dir / "*.parquet").sink_parquet(parquet_out)
    print(f"Wrote alignment stats to {parquet_out}")


if __name__ == "__main__":
    auto_cli(main)
