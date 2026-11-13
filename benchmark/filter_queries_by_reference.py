"""Keep only the query reads that align to a reference genome panel.

Raw reads sampled from a viral run include host and off-target reads. This
filter aligns every query against each reference of a panel and keeps a read
when its best alignment (most exact matches over all references and both
strands) has

    block identity  = exact matches / alignment block length  > min_identity
    query coverage  = aligned read bases / read length        >= min_coverage

Each reference is doubled head-to-tail before alignment, so reads spanning
the origin of a circular genome align end to end. minimap2 runs with the
same options as get_query_alignments_minimap.py, so the filter and the
relevance alignments downstream see reads the same way.

The panel is a two-column file (accession, label), e.g.
accessions/hbv_reference_genomes.txt. Reference fastas are fetched from NCBI
into refs_dir as genotype<label>.fa when they are not already there.

Outputs, next to output_fasta:
    <output_fasta>              the passing reads, in input order, unchanged
    <stem>.reads.tsv            one row per input read: best reference,
                                identity, coverage, pass
    <stem>.per_run.tsv          pass fraction and best-reference counts per run
    <output_fasta>.meta.json    parameters and totals

Usage:
    uv run python benchmark/filter_queries_by_reference.py \
        --queries_fa constructed/sra65viral-v2/work/queries.fa \
        --references accessions/hbv_reference_genomes.txt \
        --refs_dir constructed/sra65viral-v2/work/refs \
        --output_fasta constructed/sra65viral-v2/work/queries_hbv.fa
"""

import csv
import json
import re
import subprocess
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from jsonargparse import auto_cli
from jsonargparse.typing import Path_fr

CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch_reference(accession: str, label: str, path: Path) -> None:
    """Download one reference as fasta; a file already on disk is reused."""
    if path.exists():
        return
    query = urllib.parse.urlencode(
        {"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text"}
    )
    with urllib.request.urlopen(f"{EFETCH}?{query}", timeout=60) as response:
        text = response.read().decode()
    lines = text.strip().splitlines()
    assert lines and lines[0].startswith(">"), f"no fasta returned for {accession}"
    lines[0] = f">genotype{label}_{accession}"
    tmp = path.with_suffix(".fa.tmp")
    tmp.write_text("\n".join(lines) + "\n")
    tmp.rename(path)
    time.sleep(1)  # stay under NCBI's unauthenticated rate limit


def write_doubled(ref: Path, doubled: Path) -> None:
    """Write the reference concatenated with itself (circular genome)."""
    records = list(SeqIO.parse(ref, "fasta"))
    assert len(records) == 1, f"{ref} must hold exactly one sequence"
    seq = str(records[0].seq)
    doubled.write_text(f">{records[0].id}_x2\n{seq}{seq}\n")


def align(queries_fa: Path, ref: Path, sam: Path) -> None:
    # Same options as get_query_alignments_minimap.run_minimap
    cmd = [
        "minimap2", "-a", "--eqx", "--sam-hit-only", "-t", "1",
        str(ref), str(queries_fa),
    ]  # fmt: skip
    with open(sam, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)


def best_hits(sam: Path) -> dict[str, tuple[int, int, int, str]]:
    """Per read: (exact matches, block length, aligned read bases, strand)
    of its alignment with the most exact matches."""
    best: dict[str, tuple[int, int, int, str]] = {}
    with open(sam) as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.split("\t", 6)
            flag = int(fields[1])
            if flag & 0x4:
                continue
            eq = block = aligned = 0
            for length, op in CIGAR_OP.findall(fields[5]):
                n = int(length)
                if op == "=":
                    eq += n
                    block += n
                    aligned += n
                elif op in "XI":  # consume the read
                    block += n
                    aligned += n
                elif op in "MDN":
                    block += n
                    if op == "M":
                        aligned += n
            current = best.get(fields[0])
            if current is None or eq > current[0]:
                best[fields[0]] = (eq, block, aligned, "-" if flag & 0x10 else "+")
    return best


def main(
    queries_fa: Path_fr,
    references: Path_fr,
    refs_dir: str,
    output_fasta: str,
    min_identity: float = 0.9,
    min_coverage: float = 0.9,
):
    """
    Args:
        queries_fa: reads to filter; ids are <run accession>.<spot>.
        references: two-column file of reference accession and label.
        refs_dir: holds genotype<label>.fa (fetched when missing), the
            doubled references and the SAMs.
        output_fasta: the passing reads.
        min_identity: best-hit block identity must exceed this.
        min_coverage: best hit must cover at least this fraction of the read.
    """
    queries_fa = Path(queries_fa).resolve()
    output_fasta = Path(output_fasta).resolve()
    refs_dir = Path(refs_dir).resolve()
    refs_dir.mkdir(parents=True, exist_ok=True)
    panel = [line.split() for line in Path(references).read_text().splitlines() if line.strip()]

    # label -> best alignment per read against that reference
    hits = {}
    for accession, label in panel:
        ref = refs_dir / f"genotype{label}.fa"
        fetch_reference(accession, label, ref)
        doubled = refs_dir / f"genotype{label}.x2.fa"
        write_doubled(ref, doubled)
        sam = refs_dir / f"genotype{label}.sam"
        align(queries_fa, doubled, sam)
        hits[label] = best_hits(sam)

    rows = []
    passing = set()
    for record in SeqIO.parse(queries_fa, "fasta"):
        read_len = len(record.seq)
        row = {
            "read_id": record.id,
            "run_accession": record.id.split(".")[0],
            "read_length": read_len,
            "best_reference": "",
            "strand": "",
            "identity": 0.0,
            "coverage": 0.0,
            "matches_over_read_length": 0.0,
            "pass": False,
        }
        candidates = [
            (hits[label][record.id], label)
            for _, label in panel
            if record.id in hits[label]
        ]
        if candidates:
            # most exact matches; ties go to the earlier panel entry
            (eq, block, aligned, strand), label = max(candidates, key=lambda c: c[0][0])
            row.update(
                best_reference=label,
                strand=strand,
                identity=round(eq / block, 5),
                coverage=round(aligned / read_len, 5),
                matches_over_read_length=round(eq / read_len, 5),
            )
            row["pass"] = eq / block > min_identity and aligned / read_len >= min_coverage
        if row["pass"]:
            passing.add(record.id)
        rows.append(row)

    # Copy passing records verbatim so the fasta stays a strict subset
    keep = False
    with open(queries_fa) as src, open(output_fasta, "w") as out:
        for line in src:
            if line.startswith(">"):
                keep = line[1:].split()[0] in passing
            if keep:
                out.write(line)

    stem = output_fasta.with_suffix("")
    with open(f"{stem}.reads.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    labels = [label for _, label in panel]
    per_run = []
    for run in sorted({r["run_accession"] for r in rows}):
        run_rows = [r for r in rows if r["run_accession"] == run]
        kept = [r for r in run_rows if r["pass"]]
        best = Counter(r["best_reference"] for r in kept)
        per_run.append(
            {
                "run_accession": run,
                "reads": len(run_rows),
                "aligned": sum(1 for r in run_rows if r["best_reference"]),
                "passed": len(kept),
                "pass_fraction": round(len(kept) / len(run_rows), 4),
                **{f"passed_best_{label}": best[label] for label in labels},
            }
        )
    with open(f"{stem}.per_run.tsv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(per_run[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(per_run)

    meta = {
        "script": "benchmark/filter_queries_by_reference.py",
        "queries_fa": str(queries_fa),
        "references": {label: accession for accession, label in panel},
        "min_identity": min_identity,
        "min_coverage": min_coverage,
        "num_input": len(rows),
        "num_passed": len(passing),
    }
    Path(f"{output_fasta}.meta.json").write_text(json.dumps(meta, indent=2) + "\n")

    for r in per_run:
        print("\t".join(str(v) for v in r.values()))
    print(f"Kept {len(passing)}/{len(rows)} reads -> {output_fasta}")


if __name__ == "__main__":
    auto_cli(main, as_positional=False)
