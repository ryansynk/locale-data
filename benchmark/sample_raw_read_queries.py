import random
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO
from jsonargparse import auto_cli
from tqdm import tqdm


def stream_reads(srr_id: str, topdir: Path, num_reads: int, max_attempts: int = 4):
    """Stream the first num_reads spots of an SRA run with fastq-dump.

    fastq-dump -X streams over HTTPS and stops after num_reads spots, so the
    full run is never downloaded. NCBI's resolver service fails transiently
    ("Failed to call external services"), so retry with backoff. Accessions
    that already have a fastq on disk are skipped, making reruns resumable.
    """
    target_dir = Path(topdir) / srr_id
    target_dir.mkdir(parents=True, exist_ok=True)
    if any(target_dir.glob("*.fastq")):
        return
    for attempt in range(1, max_attempts + 1):
        try:
            subprocess.run(
                ["fastq-dump", "-X", str(num_reads), "--split-files", srr_id],
                cwd=target_dir,
                check=True,
                capture_output=True,
                text=True,
            )
            return
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr if e.stderr else str(e)
            print(f"{srr_id} attempt {attempt}/{max_attempts} failed: {error_msg}")
            if attempt == max_attempts:
                raise
            time.sleep(min(10 * 2 ** (attempt - 1), 60))
        except FileNotFoundError:
            print("SRA Toolkit not found. Is fastq-dump in your PATH?")
            raise


def get_fastq_file(srr_dir: Path) -> Path:
    """Pick the fastq to sample from, preferring mate 1 for paired runs."""
    srr_id = srr_dir.stem
    for name in (f"{srr_id}_1.fastq", f"{srr_id}.fastq", f"{srr_id}_2.fastq"):
        if (srr_dir / name).exists():
            return srr_dir / name
    raise ValueError(f"No fastq found for id {srr_id}")


def main(
    accessions_list: str,
    output_dir: str,
    output_fasta: str,
    num_reads: int = 1000,
    queries_per_accession: int | None = None,
    num_download_procs: int = 4,
    seed: int = 0,
):
    """Stream reads from SRA for a list of accessions and build an initial query set.

    Streams the first num_reads spots of each run into per-accession dirs under
    output_dir, then samples queries_per_accession reads from each (all streamed
    reads when unset) and writes them to a single query fasta. Read ids keep the
    fastq-dump format (<accession>.<spot>), which downstream scripts rely on to
    recover the source accession. Only the listed accessions are sampled
    (output_dir may hold runs from earlier sets), in sorted order, so the
    fasta is a function of the list, num_reads, and seed.
    """
    accessions_list: Path = Path(accessions_list).resolve()
    output_dir: Path = Path(output_dir).resolve()
    output_fasta: Path = Path(output_fasta).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(accessions_list) as f:
        accs = f.read().splitlines()

    # Use fewer cores for downloading to not saturate bandwidth
    with ProcessPoolExecutor(max_workers=num_download_procs) as executor:
        futures = [
            executor.submit(stream_reads, srr_id, output_dir, num_reads)
            for srr_id in accs
        ]

        for future in tqdm(
            as_completed(futures),
            desc="Streaming raw reads...",
            total=len(futures),
        ):
            future.result()

    rng = random.Random(seed)
    records = []
    for acc in tqdm(sorted(set(accs)), desc="Sampling queries..."):
        acc_records = list(SeqIO.parse(get_fastq_file(output_dir / acc), "fastq"))
        if queries_per_accession is not None and queries_per_accession < len(
            acc_records
        ):
            acc_records = rng.sample(acc_records, queries_per_accession)
        records.extend(acc_records)

    SeqIO.write(records, output_fasta, "fasta")
    print(f"Wrote {len(records)} records to {output_fasta}")


if __name__ == "__main__":
    auto_cli(main)
