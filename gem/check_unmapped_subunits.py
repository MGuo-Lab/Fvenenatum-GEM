"""Align every unmapped complex subunit (OLD_<template gene>) of the Curated model to the F. venenatum proteome.

A full Smith-Waterman scan (BLOSUM62, gap open 11, extend 1) with no prefilter, reported against the
thresholds of the orthology step (identity >= 40%, query coverage >= 70%). build_release_model.py uses the
result to decide whether a subunit is replaced or dropped.

Needs biopython. Output: model/curation/unmapped_subunit_hits.csv

Usage:  python check_unmapped_subunits.py
"""

from __future__ import annotations

import csv
import os
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from Bio.Align import PairwiseAligner, substitution_matrices

HERE = Path(__file__).resolve().parent
MODEL = HERE / "build" / "curated.xml"                  # written by curate_model.py
TARGET = HERE / "raw_data" / "fven_proteins.faa"
QUERY = HERE / "raw_data" / "aniger_cbs513.faa"
OUT = HERE.parent / "model" / "curation" / "unmapped_subunit_hits.csv"
MIN_IDENTITY, MIN_COVERAGE = 40.0, 70.0
TOP = 3
NON_RESIDUE = re.compile(r"[^ACDEFGHIKLMNPQRSTVWY]")


def read_fasta(path: Path):
    seqs, heads, key = {}, {}, None
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                head = line[1:].strip()
                m = re.match(r"\w\w\|([^|]+)\|", head)
                key = m.group(1) if m else head.split()[0]
                heads[key], seqs[key] = head, []
            elif key:
                seqs[key].append(line.strip())
    return {k: NON_RESIDUE.sub("A", "".join(v)) for k, v in seqs.items()}, heads


TARGETS, HEADS = read_fasta(TARGET)
QUERIES, _ = read_fasta(QUERY)


def aligner() -> PairwiseAligner:
    return PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                           open_gap_score=-11, extend_gap_score=-1)


def stats(al: PairwiseAligner, query: str, target: str) -> tuple[float, float]:
    """Percent identity over the alignment, gaps included, and percent coverage of the query."""
    aln = al.align(query, target)[0]
    q_blocks, t_blocks = aln.aligned
    matches = sum(sum(a == b for a, b in zip(query[q0:q1], target[t0:t1]))
                  for (q0, q1), (t0, t1) in zip(q_blocks, t_blocks))
    aligned = sum(q1 - q0 for q0, q1 in q_blocks)
    gaps = (q_blocks[-1][1] - q_blocks[0][0] - aligned) + (t_blocks[-1][1] - t_blocks[0][0] - aligned)
    return 100 * matches / (aligned + gaps), 100 * (q_blocks[-1][1] - q_blocks[0][0]) / len(query)


def scan(gene: str) -> list[dict]:
    al = aligner()
    query = QUERIES[gene]
    scored = sorted(((al.score(query, t), k) for k, t in TARGETS.items() if t), reverse=True)[:TOP]
    rows = []
    for rank, (score, hit) in enumerate(scored, 1):
        identity, coverage = stats(al, query, TARGETS[hit])
        rows.append({"template_gene": gene, "query_aa": len(query), "rank": rank, "hit": hit,
                     "hit_name": re.sub(r"^\S+ ", "", HEADS[hit]).split(" OS=")[0],
                     "identity_pct": round(identity, 1), "query_coverage_pct": round(coverage, 1),
                     "score": round(score), "meets_thresholds": identity >= MIN_IDENTITY and coverage >= MIN_COVERAGE})
    return rows


def main() -> None:
    sbml = MODEL.read_text(encoding="utf-8")
    genes = sorted(set(re.findall(r'fbc:id="G_OLD_([^"]+)"', sbml)))
    with ProcessPoolExecutor(max_workers=min(len(genes), os.cpu_count() or 2)) as pool:
        results = list(pool.map(scan, genes))
    rows = [r for group in results for r in group]
    with open(OUT, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(genes)} subunits scanned against {len(TARGETS)} proteins: {OUT}")


if __name__ == "__main__":
    main()
