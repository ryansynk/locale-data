"""Build a negative-query bundle: raw reads from runs outside an indexed set
that do not align to it, in the layout the benchmark consumes.

Searching these against the indexed set gives a null distribution of scores
(no relevant accession exists), from which false-positive-rate cutoffs are
calibrated for each method.

Reuses a superset build's work directory (e.g. constructed/sra4571/work):
its streamed raw reads (queries.fa, sample_raw_read_queries.py) and the
per-accession best-alignment parquets (best_alignments/<acc>.parquet, from
get_query_alignments_minimap.py). Each contig file is aligned on its own, so
the parquets of the indexed accessions are exactly what aligning the reads
against the indexed set alone would produce; no new alignment is needed.

Selection, in order:
  1. pool runs: in pool_accessions but not in the indexed set, whose organism
     and taxid both occur in no indexed run, with platform in `platforms` and
     assay type not in `exclude_assay_types`;
  2. per read, identity to the indexed set = max over indexed accessions of
     exact matches / read length (0 when minimap2 reports no alignment);
     keep reads with identity < identity_max, and (require_no_alignment,
     the default) only reads minimap2 did not align to any indexed
     accession at all;
  3. sample final_query_num reads, weighted so the read-length histogram
     (bins of length_bin_width bp) follows the reference bundle's queries,
     with at most max_per_accession reads per run and max_per_study per SRA
     study so one large study cannot dominate the null distribution.

Writes to output_path: queries.parquet (the positive bundle's schema with
empty relevant sets and strand "+"), queries.fa, accs.txt copied from the
reference bundle, logan_accessions/<acc> symlinks to the reference bundle's
contigs, negatives_meta.parquet (per-query provenance: source run metadata,
read length, max identity and best-hit accession in the indexed set) and
build.meta.json (parameters and counts). Run mutate_queries.py afterwards
for the per-rate query files.

Usage:
    uv run python benchmark/build_negative_queries.py \\
        --reference_bundle constructed/sra500/bundle \\
        --work_dir constructed/sra4571/work \\
        --pool_accessions accessions/sra4571.txt \\
        --metadata_tsv accessions/TableS12_SRA_Public_100studies.tsv \\
        --output_path constructed/sra500_negatives/bundle
"""

import json
import os
import shutil
from pathlib import Path

import numpy as np
import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from jsonargparse import auto_cli
from jsonargparse.typing import Path_dc, Path_drw, Path_fr
from tqdm import tqdm

QUERIES_SCHEMA = {
    "accession": pl.Utf8,
    "query_id": pl.Utf8,
    "query_sequence": pl.Utf8,
    "contig_id": pl.List(pl.Utf8),
    "contig_accession": pl.List(pl.Utf8),
    "identity": pl.List(pl.Float64),
    "ratio": pl.List(pl.Float64),
    "contig_len": pl.List(pl.Int64),
    "aln_interval_contig": pl.List(
        pl.Struct(
            {"aln_start_index_contig": pl.Int64, "aln_end_index_contig": pl.Int64}
        )
    ),
    "strand": pl.Utf8,
}

def read_pool_reads(queries_fa: Path, accessions: set[str]) -> pl.DataFrame:
    """Reads of queries_fa whose source run (id prefix before '.') is in accessions.

    Minimal fasta parser: biopython is several times slower on a multi-million
    record file, and only the id and sequence are needed.
    """
    ids, accs, seqs = [], [], []
    keep = False
    chunks: list[str] = []
    with open(queries_fa) as f:
        for line in tqdm(f, desc="Reading raw reads", unit="lines", unit_scale=True):
            if line.startswith(">"):
                if keep:
                    seqs.append("".join(chunks))
                read_id = line[1:].split(None, 1)[0]
                acc = read_id.split(".")[0]
                keep = acc in accessions
                if keep:
                    ids.append(read_id)
                    accs.append(acc)
                    chunks = []
            elif keep:
                chunks.append(line.strip())
        if keep:
            seqs.append("".join(chunks))
    return pl.DataFrame(
        {"read_id": ids, "accession": accs, "query_sequence": seqs},
        schema={"read_id": pl.Utf8, "accession": pl.Utf8, "query_sequence": pl.Utf8},
    )


def max_identity_to_indexed(
    best_alignments_dir: Path, indexed: list[str], reads: pl.DataFrame
) -> pl.DataFrame:
    """Per read: best exact-match count over the indexed accessions.

    Columns: read_id, max_matches, best_indexed_accession, n_indexed_hits.
    Reads without any alignment to an indexed accession are absent.
    """
    paths = [best_alignments_dir / f"{acc}.parquet" for acc in indexed]
    if missing := [p.name for p in paths if not p.exists()]:
        raise FileNotFoundError(
            f"{len(missing)} indexed accessions have no best-alignment parquet "
            f"under {best_alignments_dir} (e.g. {missing[:3]}); the work directory "
            "must come from a build whose accession set contains the indexed set"
        )
    # Streamed so the per-file frames are never all resident at once
    return (
        pl.scan_parquet(paths)
        .select("read_id", "accession", "matches")
        .filter(pl.col("read_id").is_in(reads["read_id"].implode()))
        .group_by("read_id")
        .agg(
            pl.col("matches").max().alias("max_matches"),
            pl.col("accession").sort_by("matches", "accession").last().alias("best_indexed_accession"),
            pl.len().alias("n_indexed_hits"),
        )
        .collect(engine="streaming")
    )


def sample_length_matched(
    candidates: pl.DataFrame,
    reference_lengths: np.ndarray,
    n: int,
    bin_width: int,
    max_per_accession: int,
    max_per_study: int,
    seed: int,
) -> pl.DataFrame:
    """Weighted sample without replacement whose read-length histogram follows
    reference_lengths, honouring per-run and per-study caps.

    Each candidate's weight is (reference mass of its length bin) / (number of
    candidates in that bin), so bins are drawn in the reference's proportions;
    a bin with no reference mass is never drawn. Sampling uses the
    Efraimidis-Spirakis keys u^(1/w): candidates are visited in key order and
    accepted unless their run or study cap is full.
    """
    ref_bins = reference_lengths // bin_width
    ref_mass = {
        int(b): c / len(ref_bins)
        for b, c in zip(*np.unique(ref_bins, return_counts=True))
    }
    cand_bins = (candidates["read_length"] // bin_width).to_numpy()
    cand_count = dict(zip(*np.unique(cand_bins, return_counts=True)))
    weights = np.array(
        [ref_mass.get(int(b), 0.0) / cand_count[b] for b in cand_bins], dtype=float
    )
    rng = np.random.default_rng(seed)
    keys = np.full(len(weights), -1.0)
    positive = weights > 0
    keys[positive] = rng.random(positive.sum()) ** (1.0 / weights[positive])

    per_acc: dict[str, int] = {}
    per_study: dict[str, int] = {}
    chosen = []
    accs = candidates["accession"].to_list()
    studies = candidates["sra_study"].to_list()
    for i in np.argsort(-keys, kind="stable"):
        if keys[i] < 0 or len(chosen) >= n:
            break
        acc, study = accs[i], studies[i]
        if (
            per_acc.get(acc, 0) >= max_per_accession
            or per_study.get(study, 0) >= max_per_study
        ):
            continue
        per_acc[acc] = per_acc.get(acc, 0) + 1
        per_study[study] = per_study.get(study, 0) + 1
        chosen.append(int(i))
    if len(chosen) < n:
        raise SystemExit(
            f"Only {len(chosen)} negatives satisfy the caps and length bins; "
            f"requested {n}. Relax max_per_accession/max_per_study or widen the pool."
        )
    return candidates[chosen]


def length_histogram(lengths: np.ndarray, bin_width: int) -> dict[str, float]:
    bins, counts = np.unique(lengths // bin_width, return_counts=True)
    return {
        f"{int(b) * bin_width}-{(int(b) + 1) * bin_width - 1}": round(
            c / len(lengths), 3
        )
        for b, c in zip(bins, counts)
    }


def link_contigs(reference_bundle: Path, output_path: Path, accessions: list[str]):
    """logan_accessions/<acc>/<acc>.contigs.fa -> the reference bundle's contig
    file, so run_benchmark finds the same contigs and never re-downloads.
    Each <acc> is a real directory holding a file symlink: Path.rglob does not
    descend into symlinked directories."""
    out_dir = output_path / "logan_accessions"
    out_dir.mkdir(exist_ok=True)
    for acc in accessions:
        target = (reference_bundle / "logan_accessions" / acc / f"{acc}.contigs.fa").resolve()
        if not target.exists():
            raise FileNotFoundError(f"reference bundle has no contigs for {acc}: {target}")
        acc_dir = out_dir / acc
        acc_dir.mkdir(exist_ok=True)
        link = acc_dir / target.name
        if link.is_symlink() or link.exists():
            link.unlink()
        os.symlink(target, link)


def main(
    reference_bundle: Path_drw,
    work_dir: Path_drw,
    pool_accessions: Path_fr,
    metadata_tsv: Path_fr,
    output_path: Path_dc,
    identity_max: float = 0.9,
    require_no_alignment: bool = True,
    final_query_num: int = 500,
    platforms: list[str] = ["ILLUMINA"],
    exclude_assay_types: list[str] = ["AMPLICON"],
    max_per_accession: int = 5,
    max_per_study: int = 50,
    length_bin_width: int = 25,
    seed: int = 0,
):
    """
    Args:
        reference_bundle: the positive bundle whose accessions are indexed
            (accs.txt, logan_accessions/, queries.parquet for the length
            distribution), e.g. constructed/sra500/bundle.
        work_dir: work directory of a build whose accession set contains the
            reference's, holding queries.fa and best_alignments/, e.g.
            constructed/sra4571/work.
        pool_accessions: accession list of that superset build; negatives are
            drawn from its runs outside the reference set.
        metadata_tsv: SRA metadata table with acc, organism, taxid, sra_study,
            platform and assay_type columns.
        output_path: bundle directory to write.
        identity_max: keep reads whose best identity (exact matches / read
            length) to any indexed accession is below this.
        require_no_alignment: additionally keep only reads with no minimap2
            alignment to any indexed accession, so no near-miss (e.g. 80%
            identity) read enters the null distribution.
        final_query_num: number of negatives to write.
        platforms: sequencing platforms admitted to the pool.
        exclude_assay_types: assay types removed from the pool.
        max_per_accession: cap on negatives drawn from one run.
        max_per_study: cap on negatives drawn from one SRA study.
        length_bin_width: bin width (bp) of the read-length histogram matched
            to the reference queries.
        seed: sampling seed.
    """
    reference_bundle = Path(reference_bundle).resolve()
    work_dir = Path(work_dir).resolve()
    output_path = Path(output_path).resolve()
    pool_accessions = Path(pool_accessions).resolve()
    metadata_tsv = Path(metadata_tsv).resolve()

    indexed = (reference_bundle / "accs.txt").read_text().split()
    pool = pool_accessions.read_text().split()
    assert len(indexed) == len(set(indexed)), "duplicate indexed accessions"
    indexed_set = set(indexed)
    outside = sorted(set(pool) - indexed_set)
    if len(outside) == len(set(pool)):
        raise SystemExit(
            "pool_accessions shares no run with the reference bundle; "
            "the work directory must have aligned against the indexed set"
        )

    # 1. Pool runs: outside the indexed set, organism and taxid unseen there.
    meta = pl.read_csv(
        metadata_tsv,
        separator="\t",
        infer_schema_length=100_000,
        columns=["acc", "organism", "taxid", "sra_study", "platform", "assay_type"],
    ).with_columns(pl.col("taxid").cast(pl.Float64).cast(pl.Int64, strict=False))
    indexed_meta = meta.filter(pl.col("acc").is_in(indexed))
    if missing := indexed_set - set(indexed_meta["acc"]):
        raise SystemExit(
            f"{len(missing)} indexed accessions missing from {metadata_tsv}"
        )
    indexed_organisms = set(indexed_meta["organism"].drop_nulls())
    indexed_taxids = set(indexed_meta["taxid"].drop_nulls())
    outside_meta = meta.filter(pl.col("acc").is_in(outside))
    n_outside = outside_meta.height
    unseen_meta = outside_meta.filter(
        ~pl.col("organism").is_in(list(indexed_organisms))
        & ~pl.col("taxid").is_in(list(indexed_taxids))
    )
    n_unseen = unseen_meta.height
    pool_meta = unseen_meta.filter(
        pl.col("platform").is_in(platforms)
        & ~pl.col("assay_type").is_in(exclude_assay_types)
    )
    pool_runs = set(pool_meta["acc"])
    print(
        f"Pool: {n_outside} runs outside the indexed set, {n_unseen} with unseen "
        f"organism+taxid, {len(pool_runs)} after platform/assay filter "
        f"({pool_meta['organism'].n_unique()} organisms, "
        f"{pool_meta['sra_study'].n_unique()} studies)"
    )

    # 2. Reads of the pool runs and their best identity to the indexed set.
    reads = read_pool_reads(work_dir / "queries.fa", pool_runs).with_columns(
        pl.col("query_sequence").str.len_chars().alias("read_length")
    )
    print(f"{reads.height} raw reads from {reads['accession'].n_unique()} pool runs")
    hits = max_identity_to_indexed(work_dir / "best_alignments", indexed, reads)
    reads = (
        reads.join(hits, on="read_id", how="left")
        .with_columns(
            (pl.col("max_matches").fill_null(0) / pl.col("read_length")).alias(
                "max_identity"
            ),
            pl.col("n_indexed_hits").fill_null(0),
        )
        .join(pool_meta.rename({"acc": "accession"}), on="accession")
    )
    candidates = reads.filter(pl.col("max_identity") < identity_max)
    print(
        f"{candidates.height} reads with identity < {identity_max} to every indexed "
        f"accession ({reads.height - candidates.height} dropped)"
    )
    if require_no_alignment:
        n_before = candidates.height
        candidates = candidates.filter(pl.col("n_indexed_hits") == 0)
        print(
            f"{candidates.height} of those with no alignment to any indexed "
            f"accession ({n_before - candidates.height} near-miss reads dropped)"
        )

    # 3. Length-matched, capped sample.
    reference_lengths = (
        pl.read_parquet(
            reference_bundle / "queries.parquet", columns=["query_sequence"]
        )
        .select(pl.col("query_sequence").str.len_chars())
        .to_series()
        .to_numpy()
    )
    chosen = sample_length_matched(
        candidates,
        reference_lengths,
        final_query_num,
        length_bin_width,
        max_per_accession,
        max_per_study,
        seed,
    )
    # Same final ordering convention as finalize_query_dataset.py
    chosen = chosen.sort("read_id").sample(len(chosen), shuffle=True, seed=seed)

    # Outputs
    output_path.mkdir(parents=True, exist_ok=True)
    queries = pl.DataFrame(
        {
            "accession": chosen["accession"],
            "query_id": chosen["read_id"],
            "query_sequence": chosen["query_sequence"],
            "contig_id": [[]] * len(chosen),
            "contig_accession": [[]] * len(chosen),
            "identity": [[]] * len(chosen),
            "ratio": [[]] * len(chosen),
            "contig_len": [[]] * len(chosen),
            "aln_interval_contig": [[]] * len(chosen),
            "strand": ["+"] * len(chosen),
        },
        schema=QUERIES_SCHEMA,
    )
    queries.write_parquet(output_path / "queries.parquet")
    SeqIO.write(
        (
            SeqRecord(Seq(seq), id=qid, description="")
            for qid, seq in queries.select("query_id", "query_sequence").iter_rows()
        ),
        output_path / "queries.fa",
        "fasta",
    )
    shutil.copy(reference_bundle / "accs.txt", output_path / "accs.txt")
    link_contigs(reference_bundle, output_path, indexed)
    chosen.select(
        pl.col("read_id").alias("query_id"),
        "accession",
        "sra_study",
        "organism",
        "taxid",
        "assay_type",
        "platform",
        "read_length",
        "max_identity",
        "best_indexed_accession",
        "n_indexed_hits",
    ).write_parquet(output_path / "negatives_meta.parquet")

    # Reads minimap2 never aligned, clear misses, and near-misses under the cutoff
    mi = chosen["max_identity"]
    strata = {
        "no alignment": int((mi == 0.0).sum()),
        "<0.8": int(((mi > 0.0) & (mi < 0.8)).sum()),
        ">=0.8": int((mi >= 0.8).sum()),
    }
    summary = {
        "script": "benchmark/build_negative_queries.py",
        "reference_bundle": str(reference_bundle),
        "work_dir": str(work_dir),
        "pool_accessions": str(pool_accessions),
        "metadata_tsv": str(metadata_tsv),
        "identity_max": identity_max,
        "require_no_alignment": require_no_alignment,
        "final_query_num": final_query_num,
        "platforms": platforms,
        "exclude_assay_types": exclude_assay_types,
        "max_per_accession": max_per_accession,
        "max_per_study": max_per_study,
        "length_bin_width": length_bin_width,
        "seed": seed,
        "counts": {
            "runs_outside_indexed": n_outside,
            "runs_unseen_organism": n_unseen,
            "runs_in_pool": len(pool_runs),
            "pool_reads": reads.height,
            "reads_eligible": candidates.height,
            "chosen": len(chosen),
            "chosen_runs": chosen["accession"].n_unique(),
            "chosen_studies": chosen["sra_study"].n_unique(),
            "chosen_organisms": chosen["organism"].n_unique(),
            "chosen_identity_strata": strata,
        },
        "length_histogram": {
            "reference": length_histogram(reference_lengths, length_bin_width),
            "chosen": length_histogram(
                chosen["read_length"].to_numpy(), length_bin_width
            ),
        },
    }
    (output_path / "build.meta.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["counts"], indent=2))
    print("length histogram (reference vs chosen):")
    for key in sorted(
        set(summary["length_histogram"]["reference"])
        | set(summary["length_histogram"]["chosen"]),
        key=lambda s: int(s.split("-")[0]),
    ):
        ref = summary["length_histogram"]["reference"].get(key, 0.0)
        got = summary["length_histogram"]["chosen"].get(key, 0.0)
        print(f"  {key:>9} bp  {ref:.3f}  {got:.3f}")
    print(
        chosen.group_by("organism")
        .agg(pl.len().alias("n"), pl.col("sra_study").n_unique().alias("studies"))
        .sort("n", descending=True)
    )
    print(f"Wrote {len(chosen)} negatives to {output_path}")


if __name__ == "__main__":
    auto_cli(main, as_positional=False)
