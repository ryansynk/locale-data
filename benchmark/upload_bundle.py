"""Upload a benchmark bundle directory to a Hugging Face dataset repo, as-is.

Every file under the bundle directory is uploaded at the same relative path
(accs.txt, queries.parquet, queries_mut<rate>.parquet, queries.fa, mutations/,
...). The dataset card is a separate tracked file (cards/<set>.md) uploaded
as README.md. All files go in one commit, so the repo never shows a partial
bundle.

Usage:
    uv run python benchmark/upload_bundle.py --bundle constructed/sra50/bundle \
        --repo_id rsynk/locale-benchmark-sra50-v2 --card cards/sra50.md
    (add --dry_run=true to list what would be uploaded without touching HF)
"""

from pathlib import Path

from huggingface_hub import CommitOperationAdd, HfApi
from jsonargparse import auto_cli
from jsonargparse.typing import Path_drw, Path_fr


def main(
    bundle: Path_drw,
    repo_id: str,
    card: Path_fr,
    private: bool = False,
    commit_message: str | None = None,
    dry_run: bool = False,
):
    """
    Args:
        bundle: bundle directory (accs.txt, queries.parquet, ...).
        repo_id: HF dataset repo, e.g. rsynk/locale-benchmark-sra50-v2
            (created if it does not exist).
        card: dataset card markdown, uploaded as README.md.
        private: create the repo private (ignored if it already exists).
        commit_message: defaults to "Upload <bundle name> bundle".
        dry_run: print the file list and exit without uploading.
    """
    bundle = Path(bundle).resolve()
    # logan_accessions/ is a tree of symlinks into the contig download that the
    # benchmark runner adds to a local bundle; it is not part of the dataset.
    files = sorted(
        p
        for p in bundle.rglob("*")
        if p.is_file()
        and not p.is_symlink()
        and "logan_accessions" not in p.relative_to(bundle).parts
    )
    if (bundle / "README.md").exists():
        raise SystemExit(
            f"{bundle} contains a README.md; the card comes from --card instead"
        )
    ops = [
        CommitOperationAdd(path_in_repo=str(p.relative_to(bundle)), path_or_fileobj=str(p))
        for p in files
    ]
    ops.append(CommitOperationAdd(path_in_repo="README.md", path_or_fileobj=str(card)))
    total = sum(p.stat().st_size for p in files) + Path(card).stat().st_size
    print(f"{repo_id}: {len(ops)} files, {total / 1e6:.1f} MB")
    for op in ops:
        print(f"  {op.path_in_repo}")
    if dry_run:
        return

    api = HfApi()
    api.create_repo(repo_id, repo_type="dataset", private=private, exist_ok=True)
    info = api.create_commit(
        repo_id=repo_id,
        repo_type="dataset",
        operations=ops,
        commit_message=commit_message or f"Upload {bundle.parent.name} bundle",
    )
    print(f"Committed: {info.commit_url}")


if __name__ == "__main__":
    auto_cli(main, as_positional=False)
