import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import polars as pl
from jsonargparse import auto_cli
from tqdm import tqdm


def logan_contigs_exist(acc: str, attempts: int = 3) -> bool:
    """Check the actual logan-pub S3 object exists, since the seqstats
    metadata can list accessions that were never synced into the release."""
    url = f"https://logan-pub.s3.amazonaws.com/c/{acc}/{acc}.contigs.fa.zst"
    req = urllib.request.Request(url, method="HEAD")
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(req, timeout=30):
                return True
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):
                return False
            raise
        except urllib.error.URLError:
            if attempt == attempts - 1:
                raise
    return False


def main(
    data_table_path: str,
    logan_metadata_path: str,
    N: int,
    output_path: str,
    seed: int = 0,
    exclude_files: list[str] | None = None,
    max_nbseq: int | None = None,
    min_bp: int | None = None,
    verify_logan: bool = True,
):
    """Sample exactly N accessions, proportionally by organism.

    This is the one step that decides which accessions a dataset may contain:
    the pool is the 100-studies table restricted to accessions the logan
    seqstats table lists (with the size filters), minus everything in
    --exclude_files. Pass the training lists there when drawing a benchmark
    set, so the output is disjoint from them by construction and nothing
    downstream has to know about train/test separation.

    With --verify_logan (the default) each draw is checked against the
    logan-pub S3 object and misses are redrawn, so the output has exactly N
    accessions that download_logan_contigs.py will succeed on. Pass
    --verify_logan=false for large training draws where a few misses are fine.

    The draw is deterministic in `seed`: organism choice, the candidate within
    the organism, and the output order.
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
    for f in exclude_files or []:
        excluded_accs.update(Path(f).read_text().split())

    # Excluded accessions leave the pool before the proportions are computed;
    # everything is sorted so the draw depends only on the seed.
    df = (
        df.filter(pl.col("acc").is_in(valid_accs) & ~pl.col("acc").is_in(excluded_accs))
        .select("acc", "organism")
        .unique(subset="acc")
        .sort("organism", "acc")
    )
    organism_proportions = (
        df.group_by("organism").len().sort("organism")
    ).with_columns((pl.col("len") / len(df)).alias("proportion"))
    organisms = organism_proportions["organism"].to_list()
    proportions = np.array(organism_proportions["proportion"].to_list())
    proportions /= proportions.sum()
    rng = np.random.default_rng(seed)

    final_accs: list[str] = []
    # Draw-verify-redraw until we have exactly N accessions that actually
    # exist in logan. `tried` prevents re-checking rejected candidates.
    tried: set[str] = set()
    num_not_in_logan = 0
    max_attempts = 50 * N
    with tqdm(total=N, desc="Sampling verified accessions...") as pbar:
        for _ in range(max_attempts):
            if len(final_accs) >= N:
                break
            s = organisms[rng.choice(len(organisms), p=proportions)]
            candidates = df.filter(
                (pl.col("organism") == s) & ~pl.col("acc").is_in(tried)
            )["acc"]
            if candidates.is_empty():
                continue
            acc = candidates[int(rng.integers(len(candidates)))]
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
    assert len(set(final_accs)) == N and not set(final_accs) & excluded_accs
    if num_not_in_logan:
        print(f"Rejected {num_not_in_logan} accessions missing from logan-pub")

    print(f"Writing {len(final_accs)} acc names to {output_path}")
    output_path.write_text("".join(acc + "\n" for acc in final_accs))


if __name__ == "__main__":
    auto_cli(main)
