"""Download reference genomes with the NCBI Datasets CLI.

Three ways to say which genomes, in priority order:

* ``--assembly_accs``: a file of assembly accessions (GCF_...), one per line.
  ``accessions/reference_genomes.txt`` is the frozen list behind the published
  reference-genome training set; this is the reproducible path.
* ``--species_list``: a file of taxon names, one per line; downloads the RefSeq
  reference assembly for each taxon that has one.
* ``--sra_accs`` (+ ``--sra_metadata``): derive the taxon list from the
  organisms of those SRA runs (what the original pipeline did with the training
  list), write it next to the zip, then download as above. Metagenome
  "organisms" simply resolve to nothing on NCBI.

``datasets`` must be on PATH or given with ``--datasets_exe``
(e.g. /pscratch/sd/r/rsynk/bin/datasets). ``--extract_to`` unzips the package
there, giving ``<extract_to>/ncbi_dataset/data`` for crop_references.py.
"""

import subprocess
import zipfile
from pathlib import Path

import polars as pl
from jsonargparse import auto_cli

DEFAULT_METADATA = Path(__file__).resolve().parents[1] / (
    "accessions/TableS12_SRA_Public_100studies.tsv"
)


def download_ncbi(
    datasets_exe: str, mode: str, input_file: Path, output_zip: Path
) -> None:
    cmd = [datasets_exe, "download", "genome", mode, "--inputfile", str(input_file)]
    if mode == "taxon":
        cmd += ["--reference", "--assembly-source", "refseq"]
    cmd += ["--include", "genome", "--filename", str(output_zip)]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def main(
    output_zip: Path,
    assembly_accs: Path | None = None,
    species_list: Path | None = None,
    sra_accs: Path | None = None,
    sra_metadata: Path = DEFAULT_METADATA,
    datasets_exe: str = "datasets",
    extract_to: Path | None = None,
):
    output_zip = Path(output_zip).resolve()
    output_zip.parent.mkdir(parents=True, exist_ok=True)

    if assembly_accs is not None:
        download_ncbi(datasets_exe, "accession", Path(assembly_accs), output_zip)
    elif species_list is not None:
        download_ncbi(datasets_exe, "taxon", Path(species_list), output_zip)
    elif sra_accs is not None:
        accs = Path(sra_accs).read_text().splitlines()
        metadata = pl.read_csv(sra_metadata, ignore_errors=True, separator="\t")
        species = (
            metadata.filter(pl.col("acc").is_in(accs))["organism"]
            .drop_nulls()
            .unique()
            .sort()
            .to_list()
        )
        species_path = output_zip.parent / "species_list.txt"
        species_path.write_text("".join(f"{s}\n" for s in species))
        print(f"wrote {len(species)} taxa to {species_path}")
        download_ncbi(datasets_exe, "taxon", species_path, output_zip)
    else:
        raise ValueError("give one of --assembly_accs, --species_list, --sra_accs")

    if extract_to is not None:
        extract_to = Path(extract_to)
        extract_to.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(output_zip) as zf:
            zf.extractall(extract_to)
        print(f"extracted to {extract_to}/ncbi_dataset/data")


if __name__ == "__main__":
    auto_cli(main)
