import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import zstandard as zstd
from jsonargparse import auto_cli
from tqdm import tqdm


def download_logan(srr_id: str, topdir: Path) -> bool:
    """Download one accession's contigs. Returns False if it's not in logan.

    An accession whose decompressed contigs.fa is already present is skipped,
    so reruns (and work dirs seeded with contigs from an earlier build) only
    fetch what is missing.
    """
    target_dir = Path(topdir) / srr_id
    fa_path = target_dir / f"{srr_id}.contigs.fa"
    if fa_path.exists() and fa_path.stat().st_size > 0:
        return True
    target_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "aws",
        "s3",
        "cp",
        f"s3://logan-pub/c/{srr_id}/{srr_id}.contigs.fa.zst",
        str(target_dir),
        "--no-sign-request",
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr if e.stderr else str(e)
        if target_dir.exists():
            shutil.rmtree(target_dir)

        # Some accessions were never synced into the logan release
        if "404" in error_msg or "does not exist" in error_msg:
            return False

        print(f"Command failed: {error_msg}")
        raise
    except FileNotFoundError:
        print("AWS CLI not found. Is 'aws' in your PATH?")
        raise
    return True


def decompress_and_cleanup(zst_path: Path):
    """Worker function to decompress a .fa.zst file and remove the original."""
    # .with_suffix('') removes the final extension (.zst), leaving .fa
    fa_path = zst_path.with_suffix("")

    # Instantiate the decompressor INSIDE the worker to avoid pickling issues
    dctx = zstd.ZstdDecompressor()

    with (
        open(zst_path, "rb") as compressed_file,
        open(fa_path, "wb") as uncompressed_file,
    ):
        # copy_stream is highly optimized for memory-efficient file decompression
        dctx.copy_stream(compressed_file, uncompressed_file)

    # Delete the original .fa.zst file
    zst_path.unlink()

    return fa_path


def main(
    accessions_list: str,
    output_dir: str,
    num_download_procs: int = 4,
    strict: bool = False,
):
    """Download logan contig sets for a candidate list of accessions.

    The input list may contain accessions not (yet) in logan; this script is
    the step that certifies the list against reality. It creates one directory
    per available accession under output_dir (each containing a verified
    non-empty <accession>.contigs.fa) and writes the true accession set to
    output_dir/downloaded_accs.txt — downstream steps should consume that file,
    not the candidate list. Pass --strict=true to instead error out when any
    requested accession is unavailable.
    """
    accessions_list: Path = Path(accessions_list).resolve()
    output_dir: Path = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(accessions_list) as f:
        accs = f.read().splitlines()

    # Use fewer cores for downloading to not saturate bandwidth
    with ProcessPoolExecutor(max_workers=num_download_procs) as download_executor:
        futures = {
            download_executor.submit(download_logan, srr_id, output_dir): srr_id
            for srr_id in accs
        }

        missing = []
        for future in tqdm(
            as_completed(futures),
            desc="Downloading logan contigs...",
            total=len(futures),
        ):
            if not future.result():
                missing.append(futures[future])

    logan_paths = list(output_dir.rglob("*.fa.zst"))
    with ProcessPoolExecutor() as decompress_executor:
        futures = [
            decompress_executor.submit(decompress_and_cleanup, logan_path)
            for logan_path in logan_paths
        ]
        for f in tqdm(
            as_completed(futures), desc="Decompressing files...", total=len(futures)
        ):
            f.result()

    # Certify: every non-missing accession must have a non-empty contigs.fa
    downloaded = []
    for srr_id in accs:
        fa_path = output_dir / srr_id / f"{srr_id}.contigs.fa"
        if fa_path.exists() and fa_path.stat().st_size > 0:
            downloaded.append(srr_id)
        elif srr_id not in missing:
            missing.append(srr_id)

    downloaded_accs_file = output_dir / "downloaded_accs.txt"
    downloaded_accs_file.write_text("".join(acc + "\n" for acc in downloaded))
    print(
        f"{len(downloaded)}/{len(accs)} accessions downloaded and verified; "
        f"certified list written to {downloaded_accs_file}"
    )
    if missing:
        print(f"Missing from logan or failed verification: {missing}")
        if strict:
            raise SystemExit(
                f"{len(missing)} of {len(accs)} accessions unavailable. "
                "Re-run subsample_accessions.py to top up the list, or pass "
                "--strict=false to proceed with the certified subset."
            )


if __name__ == "__main__":
    auto_cli(main)
