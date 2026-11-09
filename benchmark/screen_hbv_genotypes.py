"""Screen SRA runs for HBV genotype label agreement, mapping rate and breadth.

Replicates the selection screen behind the SRA-viral benchmark. For each run,
the first `num_reads` spots are streamed with fastq-dump and mapped with
minimap2 against each genotype reference separately. Per mapped read the
identity is matching bases over alignment block length; the run's label is the
reference with the highest median identity. A run also reports its mapping
rate (mapped reads / reads examined, max over references) and its breadth on
the best reference (fraction of reference bases at depth >= min_depth).

The caller decides which runs pass; this script only measures. Results are
written as a tsv so the per-project label screen and the per-run screens can
be applied without re-mapping.

Usage:
    uv run python benchmark/screen_hbv_genotypes.py \
        --runs_tsv <tsv with run_accession,study_accession,genotype> \
        --refs_dir <dir of genotype<X>.fa> \
        --work_dir <scratch for fastqs and sams> \
        --output_tsv <results.tsv>
"""

import re
import statistics
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")


def stream_reads(run: str, work_dir: Path, num_reads: int) -> Path | None:
    """Stream the first num_reads spots of a run; return the fastq to map."""
    target = work_dir / "fastq" / run
    target.mkdir(parents=True, exist_ok=True)
    existing = [f for f in sorted(target.glob("*.fastq")) if not f.name.endswith("_all.fastq")]
    if not existing:
        try:
            subprocess.run(
                ["fastq-dump", "-X", str(num_reads), "--split-files", run],
                cwd=target,
                check=True,
                capture_output=True,
                text=True,
                timeout=1800,
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            print(f"{run}: fastq-dump failed: {e}", file=sys.stderr)
            return None
        existing = sorted(target.glob("*.fastq"))
    mates = [f for f in sorted(target.glob("*.fastq")) if not f.name.endswith("_all.fastq")]
    if not mates:
        return None
    if len(mates) == 1:
        return mates[0]
    # Both mates are mapped: for tiled-amplicon runs mate 1 alone covers only
    # part of the genome, which under-measures breadth. They are mapped as
    # independent single-end reads, which is all the identity/coverage
    # statistics need.
    merged = target / f"{run}_all.fastq"
    if not merged.exists() or merged.stat().st_size == 0:
        tmp = merged.with_suffix(".fastq.tmp")
        with open(tmp, "w") as out:
            for mate in mates:
                with open(mate) as handle:
                    for line in handle:
                        out.write(line)
        tmp.rename(merged)
    return merged


def count_reads(fastq: Path) -> int:
    with open(fastq) as handle:
        return sum(1 for i, _ in enumerate(handle) if i % 4 == 0)


def map_to_ref(fastq: Path, ref: Path, sam: Path) -> None:
    if sam.exists():
        return
    tmp = sam.with_suffix(".sam.tmp")
    with open(tmp, "w") as out:
        subprocess.run(
            ["minimap2", "-ax", "sr", "--eqx", "--secondary=no", "-t", "1",
             str(ref), str(fastq)],
            stdout=out, stderr=subprocess.DEVNULL, check=True,
        )
    tmp.rename(sam)


def ref_length(ref: Path) -> int:
    return sum(len(line.strip()) for line in open(ref) if not line.startswith(">"))


def parse_sam(sam: Path, reflen: int) -> tuple[list[float], int, list[int]]:
    """Return (per-read identities, mapped read count, coverage array)."""
    identities: list[float] = []
    coverage = [0] * reflen
    mapped = 0
    with open(sam) as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.split("\t", 6)
            flag = int(fields[1])
            if flag & 0x4 or flag & 0x900:
                continue
            mapped += 1
            pos = int(fields[3]) - 1
            eq = block = 0
            ref_i = pos
            for length, op in CIGAR_OP.findall(fields[5]):
                n = int(length)
                if op == "=":
                    eq += n
                    block += n
                    for i in range(ref_i, min(ref_i + n, reflen)):
                        coverage[i] += 1
                    ref_i += n
                elif op == "X":
                    block += n
                    for i in range(ref_i, min(ref_i + n, reflen)):
                        coverage[i] += 1
                    ref_i += n
                elif op in "MDN":
                    block += n
                    for i in range(ref_i, min(ref_i + n, reflen)):
                        coverage[i] += 1
                    ref_i += n
                elif op == "I":
                    block += n
            if block:
                identities.append(eq / block)
    return identities, mapped, coverage


def screen_run(run: str, refs: list[Path], work_dir: Path, num_reads: int,
               min_depth: int) -> dict | None:
    fastq = stream_reads(run, work_dir, num_reads)
    if fastq is None:
        return None
    total = count_reads(fastq)
    if total == 0:
        return None
    sam_dir = work_dir / "sam" / run
    sam_dir.mkdir(parents=True, exist_ok=True)
    row = {"run_accession": run, "reads_examined": total}
    best = (None, -1.0)
    per_ref = {}
    for ref in refs:
        geno = ref.stem.replace("genotype", "")
        sam = sam_dir / f"{geno}.sam"
        map_to_ref(fastq, ref, sam)
        reflen = ref_length(ref)
        identities, mapped, coverage = parse_sam(sam, reflen)
        med = statistics.median(identities) if identities else 0.0
        rate = mapped / total
        breadth = sum(1 for c in coverage if c >= min_depth) / reflen
        per_ref[geno] = (med, rate, breadth)
        row[f"median_identity_{geno}"] = round(med, 5)
        row[f"map_rate_{geno}"] = round(rate, 5)
        row[f"breadth_{geno}"] = round(breadth, 5)
        if med > best[1]:
            best = (geno, med)
    label = best[0]
    row["label"] = label
    row["label_median_identity"] = round(per_ref[label][0], 5)
    row["best_map_rate"] = round(max(v[1] for v in per_ref.values()), 5)
    row["label_map_rate"] = round(per_ref[label][1], 5)
    row["label_breadth"] = round(per_ref[label][2], 5)
    return row


def main(runs_tsv: str, refs_dir: str, work_dir: str, output_tsv: str,
         num_reads: int = 2000, min_depth: int = 10, num_workers: int = 4):
    import csv

    refs = sorted(Path(refs_dir).glob("genotype?.fa"))
    assert refs, f"no genotype<X>.fa in {refs_dir}"
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    meta = {r["run_accession"]: r for r in csv.DictReader(open(runs_tsv), delimiter="\t")}
    runs = list(meta)
    print(f"screening {len(runs)} runs against {[r.stem for r in refs]}")

    rows = []
    with ProcessPoolExecutor(max_workers=num_workers) as ex:
        futures = {ex.submit(screen_run, run, refs, work, num_reads, min_depth): run
                   for run in runs}
        for fut in as_completed(futures):
            run = futures[fut]
            try:
                row = fut.result()
            except Exception as e:  # noqa: BLE001
                print(f"{run}: {e}", file=sys.stderr)
                row = None
            if row is None:
                print(f"{run}: no result", file=sys.stderr)
                continue
            row["study_accession"] = meta[run].get("study_accession", "")
            row["metadata_genotype"] = meta[run].get("genotype", "")
            rows.append(row)
            print(f"  {run} label={row['label']} id={row['label_median_identity']:.3f} "
                  f"rate={row['label_map_rate']:.3f} breadth={row['label_breadth']:.3f} "
                  f"(meta {row['metadata_genotype']})", flush=True)

    if not rows:
        raise SystemExit("no runs screened")
    cols = ["run_accession", "study_accession", "metadata_genotype", "label",
            "label_median_identity", "label_map_rate", "label_breadth",
            "best_map_rate", "reads_examined"]
    cols += [c for c in rows[0] if c not in cols]
    with open(output_tsv, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        writer.writeheader()
        for row in sorted(rows, key=lambda r: (r["metadata_genotype"], r["run_accession"])):
            writer.writerow(row)
    print(f"wrote {len(rows)} rows to {output_tsv}")


if __name__ == "__main__":
    from jsonargparse import auto_cli

    auto_cli(main, as_positional=False)
