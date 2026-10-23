"""Score a result file against the LOCALE benchmark bundle.

Reports two families of numbers side by side:

  vector recall@k   -- the standard ANN metric. How many of the exact
                       nearest neighbours your index actually returned.
  accession metrics -- the downstream task. Retrieved vector ids are mapped
                       to SRA accessions and scored against the accessions the
                       query's read actually aligns to.

The gap between them is the interesting part: you only need one chunk of the
right accession in your top-k, so accession recall stays high long after
vector recall degrades.

Result file format (identical to gt.bin, so one reader serves both):

    uint32   nq
    uint32   k
    uint32   ids[nq * k]      row indices into the base embeddings
    float32  dists[nq * k]    inner-product scores, higher is better

Rows are query *chunks*, not queries. query_meta.parquet maps chunks back to
query ids; a query longer than max_seq_len owns several consecutive rows.

Depends on numpy and pyarrow only -- no torch, no genomics stack.
"""

import argparse
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

# Accessions the method never returned are left here so they can be told apart
# from a genuine low score. They must not be credited when a ranked cut lands
# inside the tied block of unscored accessions.
DEFAULT_SCORE = -2.0


# --------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------


def read_result_bin(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a gt.bin / res.bin into (ids, dists), both (nq, k)."""
    with open(path, "rb") as f:
        nq, k = np.frombuffer(f.read(8), dtype=np.uint32)
        nq, k = int(nq), int(k)
        ids = np.frombuffer(f.read(nq * k * 4), dtype=np.uint32).reshape(nq, k)
        dists = np.frombuffer(f.read(nq * k * 4), dtype=np.float32).reshape(nq, k)
    return ids.astype(np.int64), dists.astype(np.float64)


def read_meta(path: Path) -> tuple[list[str], np.ndarray]:
    """Read meta.parquet into (accession names, start_row boundaries).

    Rows [start_row, start_row + num_rows) of the base belong to that
    accession, so a row id resolves by searchsorted on the starts.
    """
    table = pq.read_table(path)
    names = table.column("srr_id").to_pylist()
    starts = np.asarray(table.column("start_row").to_pylist(), dtype=np.int64)
    counts = np.asarray(table.column("num_rows").to_pylist(), dtype=np.int64)

    order = np.argsort(starts)
    names = [names[i] for i in order]
    starts, counts = starts[order], counts[order]
    return names, starts, counts


def read_query_meta(path: Path) -> tuple[list[str], np.ndarray]:
    table = pq.read_table(path)
    query_ids = table.column("query_id").to_pylist()
    starts = np.asarray(table.column("start_row").to_pylist(), dtype=np.int64)
    counts = np.asarray(table.column("num_rows").to_pylist(), dtype=np.int64)
    return query_ids, starts, counts


def read_labels(path: Path) -> dict[str, list[str]]:
    """query_id -> accessions the read aligns to (the relevant set)."""
    table = pq.read_table(path, columns=["query_id", "contig_accession"])
    return dict(
        zip(
            table.column("query_id").to_pylist(),
            table.column("contig_accession").to_pylist(),
        )
    )


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------


def _drop_unscored(indices: np.ndarray, y_score: np.ndarray) -> np.ndarray:
    """Drop accessions the method never scored from a ranked selection.

    Unscored accessions all sit at DEFAULT_SCORE and therefore tie. argsort
    breaks ties by index, so a cut landing inside that block would hand out
    credit for accessions the method never returned.
    """
    return indices[y_score[indices] > DEFAULT_SCORE]


def recall_at_k(y_true: np.ndarray, y_score: np.ndarray, k: int) -> float:
    n_relevant = int(y_true.sum())
    if n_relevant == 0:
        return np.nan
    top = (-y_score).argsort()[:k]
    top = _drop_unscored(top, y_score)
    return int(y_true[top].sum()) / n_relevant


def r_precision(y_true: np.ndarray, y_score: np.ndarray) -> float:
    n_relevant = int(y_true.sum())
    if n_relevant == 0:
        return np.nan
    top = (-y_score).argsort()[:n_relevant]
    top = _drop_unscored(top, y_score)
    return int(y_true[top].sum()) / n_relevant


def average_precision(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Binary average precision, matching sklearn's average_precision_score."""
    n_relevant = int(y_true.sum())
    if n_relevant == 0:
        return np.nan
    order = (-y_score).argsort()
    hits = y_true[order]
    precision = np.cumsum(hits) / np.arange(1, len(hits) + 1)
    return float(np.sum(precision * hits) / n_relevant)


def vector_recall(res_ids: np.ndarray, gt_ids: np.ndarray, k: int) -> float:
    """Mean |returned top-k ∩ true top-k| / k over query chunks."""
    hits = 0
    for res_row, gt_row in zip(res_ids[:, :k], gt_ids[:, :k]):
        hits += np.intersect1d(res_row, gt_row, assume_unique=False).size
    return hits / (res_ids.shape[0] * k)


# --------------------------------------------------------------------------
# collapse: rows -> accessions -> queries
# --------------------------------------------------------------------------


def collapse_to_accessions(
    ids: np.ndarray,
    dists: np.ndarray,
    base_starts: np.ndarray,
    n_acc: int,
    chunk_to_query: np.ndarray,
    n_queries: int,
) -> np.ndarray:
    """Build the (n_queries, n_acc) score matrix.

    An accession scores as the best-scoring chunk of it that any chunk of the
    query retrieved -- max, not sum, so a large accession is not rewarded
    merely for having more vectors in the index.
    """
    acc_idx = np.searchsorted(base_starts, ids.ravel(), side="right") - 1
    query_idx = np.repeat(chunk_to_query, ids.shape[1])

    scores = np.full((n_queries, n_acc), DEFAULT_SCORE, dtype=np.float64)
    np.maximum.at(scores, (query_idx, acc_idx), dists.ravel())
    return scores


def build_truth(
    query_ids: list[str],
    labels: dict[str, list[str]],
    acc_to_idx: dict[str, int],
) -> np.ndarray:
    truth = np.zeros((len(query_ids), len(acc_to_idx)), dtype=np.int64)
    for i, query_id in enumerate(query_ids):
        for accession in labels.get(query_id, []):
            idx = acc_to_idx.get(accession)
            if idx is not None:
                truth[i, idx] = 1
    return truth


# --------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("results", type=Path, help="your res.bin")
    parser.add_argument(
        "--bundle",
        type=Path,
        required=True,
        help="directory holding meta.parquet and the mut*/ subdirectories",
    )
    parser.add_argument(
        "--mut",
        type=float,
        default=0.00,
        help="mutation rate of the query set the results came from",
    )
    parser.add_argument(
        "--queries",
        type=Path,
        default=None,
        help="queries.parquet with the labels (default: <bundle>/../../queries.parquet)",
    )
    parser.add_argument(
        "--k",
        default="1,5,10",
        help="comma-separated accession-level cutoffs (default: 1,5,10)",
    )
    parser.add_argument("--label", default=None, help="name for this run in the output")
    args = parser.parse_args()

    bundle = args.bundle
    mut_dir = bundle / f"mut{args.mut:.2f}"
    queries_path = args.queries or bundle.parent.parent / "queries.parquet"
    ks = [int(x) for x in args.k.split(",")]

    res_ids, res_dists = read_result_bin(args.results)
    gt_ids, _ = read_result_bin(mut_dir / "gt.bin")
    acc_names, base_starts, _ = read_meta(bundle / "meta.parquet")
    query_ids, query_starts, query_counts = read_query_meta(
        mut_dir / "query_meta.parquet"
    )
    labels = read_labels(queries_path)

    if res_ids.shape[0] != gt_ids.shape[0]:
        raise SystemExit(
            f"row mismatch: results have {res_ids.shape[0]} rows, "
            f"gt.bin has {gt_ids.shape[0]}. Both must be one row per query "
            f"chunk -- see query_meta.parquet."
        )

    # chunk row -> query index
    chunk_to_query = np.repeat(np.arange(len(query_ids)), query_counts)
    if chunk_to_query.size != res_ids.shape[0]:
        raise SystemExit(
            f"query_meta.parquet describes {chunk_to_query.size} chunks but the "
            f"result file has {res_ids.shape[0]} rows"
        )

    acc_to_idx = {name: i for i, name in enumerate(acc_names)}
    scores = collapse_to_accessions(
        res_ids,
        res_dists,
        base_starts,
        len(acc_names),
        chunk_to_query,
        len(query_ids),
    )
    truth = build_truth(query_ids, labels, acc_to_idx)

    scored = [i for i in range(len(query_ids)) if truth[i].sum() > 0]
    if not scored:
        raise SystemExit("no query had a relevant accession; check --queries")

    k_res = res_ids.shape[1]
    name = args.label or args.results.name

    print(f"\n{name}")
    print(f"  queries      {len(scored)} scored of {len(query_ids)}")
    print(f"  chunks       {res_ids.shape[0]}   k={k_res}   mut={args.mut:.2f}")
    print(f"  accessions   {len(acc_names)}")
    print()
    print(f"  vector recall@{k_res}   {vector_recall(res_ids, gt_ids, k_res):.4f}")
    print()
    for k in ks:
        values = [recall_at_k(truth[i], scores[i], k) for i in scored]
        print(f"  accession recall@{k:<4} {np.nanmean(values):.4f}")
    rp = [r_precision(truth[i], scores[i]) for i in scored]
    ap = [average_precision(truth[i], scores[i]) for i in scored]
    print(f"  R-precision          {np.nanmean(rp):.4f}")
    print(f"  AUPRC                {np.nanmean(ap):.4f}")
    print()


if __name__ == "__main__":
    main()
