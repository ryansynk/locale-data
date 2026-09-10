import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import polars as pl
from jsonargparse import auto_cli
from tqdm import tqdm


def logan_contigs_exist(acc: str) -> bool:
    """Check the actual logan-pub S3 object exists, since the seqstats
    metadata can list accessions that were never synced into the release."""
    url = f"https://logan-pub.s3.amazonaws.com/c/{acc}/{acc}.contigs.fa.zst"
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=30):
            return True
    except urllib.error.HTTPError as e:
        if e.code in (403, 404):
            return False
        raise


def main(
    data_table_path: str,
    logan_metadata_path: str,
    N: int,
    output_path: str,
    exclude_files: list[str] | None = None,
    max_nbseq: int | None = None,
    min_bp: int | None = None,
    verify_logan: bool = False,
):
    """Sample N candidate accessions, proportionally by organism.

    The output is a candidate list: download_logan_contigs.py certifies it
    against what actually exists in logan and emits the true list
    (downloaded_accs.txt). Pass --verify_logan=true to instead verify each
    candidate against S3 at draw time, guaranteeing exactly N verified
    accessions (or a loud failure if the pool is exhausted).
    """
    output_path: Path = Path(output_path).resolve()
    df = pl.read_csv(data_table_path, separator="\t", ignore_errors=True)
    logan_metadata = pl.read_parquet(logan_metadata_path)
    logan_metadata = logan_metadata.with_columns(
        pl.all().exclude("accession").cast(pl.Int64)
    )
    if max_nbseq:
        logan_metadata = logan_metadata.filter(
            pl.col("seqstats_contigs_nbseq") <= max_nbseq
        )
    if min_bp:
        logan_metadata = logan_metadata.filter(
            pl.col("seqstats_contigs_sumlen") >= min_bp
        )

    valid_accs = set(logan_metadata["accession"].to_list())

    excluded_accs: set[str] = set()
    if exclude_files:
        for f in exclude_files:
            excluded_accs.update(Path(f).read_text().splitlines())

    df = df.filter(pl.col("acc").is_in(valid_accs))

    organism_proportions = (
        df.group_by("organism")
        .len()
        .with_columns((pl.col("len") / len(df)).alias("proportion"))
    )
    organisms = organism_proportions["organism"].to_list()
    proportions = organism_proportions["proportion"].to_list()
    rng = np.random.default_rng()

    final_accs: list[str] = []
    # Draw-verify-redraw until we have exactly N accessions that actually
    # exist in logan. `tried` prevents re-checking rejected candidates.
    tried: set[str] = set(excluded_accs)
    num_not_in_logan = 0
    max_attempts = 50 * N
    with tqdm(total=N, desc="Sampling verified accessions...") as pbar:
        for _ in range(max_attempts):
            if len(final_accs) >= N:
                break
            s = rng.choice(organisms, p=proportions)
            candidates = df.filter(
                (pl.col("organism") == s) & ~pl.col("acc").is_in(tried)
            )
            if candidates.is_empty():
                continue
            acc = candidates.sample(1)["acc"].item()
            tried.add(acc)
            if verify_logan and not logan_contigs_exist(acc):
                num_not_in_logan += 1
                continue
            final_accs.append(acc)
            pbar.update(1)

    assert len(final_accs) == N, (
        f"Could only sample {len(final_accs)} of {N} requested accessions "
        f"({num_not_in_logan} candidates rejected as missing from logan). "
        "The candidate pool is exhausted — relax the filters or lower N."
    )
    if num_not_in_logan:
        print(f"Rejected {num_not_in_logan} accessions missing from logan-pub")

    print(f"Writing {len(final_accs)} acc names to {output_path}")
    with open(output_path, "w") as f:
        for acc in final_accs:
            f.write(acc + "\n")


if __name__ == "__main__":
    auto_cli(main)
