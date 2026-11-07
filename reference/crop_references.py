"""Crop NCBI reference genomes into a fixed-length training parquet.

Reads an NCBI Datasets data package (``data_dir/assembly_data_report.jsonl``
plus ``data_dir/<assembly accession>/*.fna``), cuts every record into
non-overlapping ``crop_len`` windows, labels each crop with the organism name
used in the SRA metadata (so the ``organism`` column joins against
``accessions/TableS12_SRA_Public_100studies.tsv``), shuffles, and writes a
parquet with columns ``sequence / name / accession / organism / crop_start /
crop_end`` -- the same schema as the original reference_genome_dataset.parquet.

Genomes can be left out with ``--exclude_organisms`` (SRA organism names, e.g.
``"[Homo sapiens,Mus musculus]"``) or ``--exclude_accessions`` (assembly
accessions, e.g. ``"[GCF_000001405.40]"``).
"""

from pathlib import Path

import polars as pl
from jsonargparse import auto_cli
from pyfaidx import Fasta
from tqdm import tqdm

DEFAULT_DATA_DIR = Path("/pscratch/sd/r/rsynk/reference_data/refs/ncbi_dataset/data")
DEFAULT_METADATA = Path(__file__).resolve().parents[1] / (
    "accessions/TableS12_SRA_Public_100studies.tsv"
)

# Maps NCBI assembly organismName -> SRA organism name for cases where NCBI has
# updated the taxonomy (synonym renames) that can't be resolved by prefix matching.
ORGANISM_SYNONYMS: dict[str, str] = {
    "Lathyrus oleraceus": "Pisum sativum",
}


def _map_assembly_name_to_sra_name(
    assembly_name: str, sra_names: set[str]
) -> str | None:
    """Return the SRA organism name that corresponds to an assembly organismName.

    Tries exact match, then ORGANISM_SYNONYMS, then prefix matching (handles
    strain-appended names like "Bacillus anthracis str. Ames" matching SRA's
    "Bacillus anthracis").
    """
    if assembly_name in sra_names:
        return assembly_name
    if (
        assembly_name in ORGANISM_SYNONYMS
        and ORGANISM_SYNONYMS[assembly_name] in sra_names
    ):
        return ORGANISM_SYNONYMS[assembly_name]
    prefix_matches = [n for n in sra_names if assembly_name.startswith(n)]
    if prefix_matches:
        return max(prefix_matches, key=len)
    return None


def get_name_map(
    metadata_path: Path,
    assembly_organism_names: list[str],
    accs_path: Path | None = None,
) -> dict[str, str]:
    """Map each assembly organismName to its SRA organism name.

    Args:
        metadata_path: SRA metadata TSV with 'acc' and 'organism' columns.
        assembly_organism_names: organismName values from the assembly report.
        accs_path: optional SRA accession list; if given, only the organisms of
            those runs are candidates (the original pipeline used the training
            list). Default: every organism in the metadata.

    Raises if any assembly organism cannot be matched.
    """
    metadata = pl.read_csv(metadata_path, ignore_errors=True, separator="\t")
    if accs_path is not None:
        accs = Path(accs_path).read_text().splitlines()
        metadata = metadata.filter(pl.col("acc").is_in(accs))
    sra_names: set[str] = set(metadata["organism"].drop_nulls().to_list())

    name_map: dict[str, str] = {}
    unmatched = []
    for assembly_name in assembly_organism_names:
        sra_name = _map_assembly_name_to_sra_name(assembly_name, sra_names)
        if sra_name is not None:
            name_map[assembly_name] = sra_name
        else:
            unmatched.append(assembly_name)

    assert len(unmatched) == 0, (
        f"{len(unmatched)} assembly organism(s) could not be matched to "
        f"any SRA organism name:\n" + "\n".join(f"  {n}" for n in unmatched)
    )
    return name_map


def count_crops(genome_files: list[Path], crop_len: int, stride: int) -> int:
    """Number of windows get_crops would produce (from the faidx lengths only)."""
    n = 0
    for f in genome_files:
        for name, rec in Fasta(str(f)).items():
            n += max(0, (len(rec) - crop_len) // stride + 1)
    return n


def stride_for_target(
    genome_files: list[Path], crop_len: int, target_rows: int
) -> int:
    """Largest stride (<= crop_len) that yields at least target_rows windows."""
    lo, hi = 1, crop_len  # count is non-increasing in stride
    if count_crops(genome_files, crop_len, hi) >= target_rows:
        return hi
    assert count_crops(genome_files, crop_len, lo) >= target_rows, (
        f"even stride 1 gives fewer than {target_rows:,} crops"
    )
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if count_crops(genome_files, crop_len, mid) >= target_rows:
            lo = mid
        else:
            hi = mid
    return lo


def get_crops(
    genome_file: Path,
    crop_len: int,
    accession: str,
    sra_organism: str,
    stride: int | None = None,
) -> pl.DataFrame:
    """Cut every record of genome_file into crop_len windows every stride bp.

    stride defaults to crop_len (non-overlapping, the original dataset). A
    trailing partial window (< crop_len) is dropped. Sequences are uppercased;
    crops containing N are kept, as in the original dataset.
    """
    stride = crop_len if stride is None else stride
    genes = Fasta(str(genome_file))
    all_names: list[str] = []
    all_seqs: list[str] = []
    all_crop_starts: list[int] = []
    all_crop_ends: list[int] = []
    for name in genes.keys():
        seq = str(genes[name]).upper()
        for start in range(0, len(seq) - crop_len + 1, stride):
            all_names.append(name)
            all_seqs.append(seq[start : start + crop_len])
            all_crop_starts.append(start)
            all_crop_ends.append(start + crop_len)

    return pl.DataFrame(
        {
            "sequence": all_seqs,
            "name": all_names,
            "accession": pl.Series([accession] * len(all_seqs), dtype=pl.String),
            "organism": pl.Series([sra_organism] * len(all_seqs), dtype=pl.String),
            "crop_start": pl.Series(all_crop_starts, dtype=pl.Int64),
            "crop_end": pl.Series(all_crop_ends, dtype=pl.Int64),
        }
    )


def main(
    output_path: Path,
    data_dir: Path = DEFAULT_DATA_DIR,
    metadata_path: Path = DEFAULT_METADATA,
    accs_path: Path | None = None,
    crop_len: int = 2048,
    stride: int | None = None,
    target_rows: int | None = None,
    exclude_organisms: list[str] = [],
    exclude_accessions: list[str] = [],
    seed: int = 42,
):
    """Build the reference-genome training parquet.

    Args:
        output_path: parquet to write.
        data_dir: NCBI Datasets data dir (holds assembly_data_report.jsonl and
            one <assembly accession>/ subdir per genome).
        metadata_path: SRA metadata TSV ('acc', 'organism') used to label crops
            with the SRA organism name.
        accs_path: optional SRA accession list restricting the organism-name
            candidates to those runs' organisms.
        crop_len: window length in bp.
        stride: distance between window starts (default crop_len: non-overlapping).
        target_rows: if given, choose the largest stride that yields at least
            this many windows, then keep exactly this many after shuffling
            (e.g. 6768336 to match the original 7-genome parquet after
            excluding genomes). Overrides stride.
        exclude_organisms: SRA organism names (after mapping) to leave out.
        exclude_accessions: assembly accessions (GCF_...) to leave out.
        seed: shuffle seed.
    """
    data_dir = Path(data_dir)
    assembly_metadata = pl.read_ndjson(data_dir / "assembly_data_report.jsonl")
    # Key by the downloaded "accession" (= directory name), not "currentAccession",
    # which can point at a newer version (Danio rerio GCF_049306965.1 -> .2).
    assembly_organisms = assembly_metadata.unnest("organism").select(
        ["accession", "organismName"]
    )
    name_map = get_name_map(
        metadata_path,
        assembly_organisms["organismName"].unique().to_list(),
        accs_path,
    )
    accession_to_sra_organism = {
        row["accession"]: name_map[row["organismName"]]
        for row in assembly_organisms.iter_rows(named=True)
    }

    # Resolve exclusions, insisting that every requested exclusion hits something
    # (a typo would otherwise silently keep the genome).
    excluded: set[str] = set()
    for org in exclude_organisms:
        hits = {a for a, o in accession_to_sra_organism.items() if o == org}
        assert hits, (
            f"--exclude_organisms {org!r} matches no genome; available: "
            f"{sorted(set(accession_to_sra_organism.values()))}"
        )
        excluded |= hits
    for acc in exclude_accessions:
        assert acc in accession_to_sra_organism, (
            f"--exclude_accessions {acc!r} not in data package; available: "
            f"{sorted(accession_to_sra_organism)}"
        )
        excluded.add(acc)

    genome_files = sorted(data_dir.glob("*/*.fna"))  # glob, not rglob: follows symlinked dirs
    assert genome_files, f"no .fna files under {data_dir}"
    kept_files = [f for f in genome_files if f.parent.name not in excluded]
    for f in genome_files:
        acc = f.parent.name
        status = "EXCLUDED" if acc in excluded else "keep"
        print(f"{status:>8}  {acc:<18} {accession_to_sra_organism[acc]}")

    if target_rows is not None:
        stride = stride_for_target(kept_files, crop_len, target_rows)
        n = count_crops(kept_files, crop_len, stride)
        print(f"stride {stride} -> {n:,} windows, keeping {target_rows:,}")
    elif stride is None:
        stride = crop_len

    all_crops = []
    for f in tqdm(kept_files, desc="Cropping genomes"):
        accession = f.parent.name
        all_crops.append(
            get_crops(
                f, crop_len, accession, accession_to_sra_organism[accession], stride
            )
        )
    crops_df = pl.concat(all_crops).sample(fraction=1.0, shuffle=True, seed=seed)
    if target_rows is not None:
        crops_df = crops_df.head(target_rows)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    crops_df.write_parquet(output_path)
    print(f"wrote {len(crops_df):,} crops of {crop_len} bp (stride {stride}) to {output_path}")
    print(crops_df.group_by("accession", "organism").len().sort("len"))


if __name__ == "__main__":
    auto_cli(main)
