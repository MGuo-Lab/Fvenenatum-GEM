"""Merge OrthoFinder pairwise orthologs with the existing RBH orthologs.

Strategy
========
- OrthoFinder gives tree-based orthology inferred across all 10 species in
  the input set, which captures distant orthologs that pure-RBH BLAST misses
  (orthologs that fall below the ≥40% identity / ≥70% coverage cutoff but
  share gene-tree topology with the F. venenatum gene).

- RBH gives high-confidence direct ortholog evidence at strict thresholds.
  Where both agree, that's a doubly-supported call.

The merged table keeps every template gene with at least one
F. venenatum ortholog, with a ``method`` column flagging whether the call
came from RBH, OrthoFinder, or both.

For multi-mapped iJB1325 genes (multiple F. venenatum candidates), pick the
one with the highest BLAST percent identity if available; otherwise the
first OrthoFinder ortholog.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
from collections import defaultdict


def load_rbh(path: str) -> dict[str, dict]:
    """Load existing RBH orthologs.tsv → {template_gene: {target_gene, pident, evalue}}"""
    rbh: dict[str, dict] = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            t = row["template_gene"].strip()
            if t:
                rbh[t] = {
                    "target":     row["target_gene"].strip(),
                    "pident_fwd": row.get("pident_fwd", ""),
                    "pident_rev": row.get("pident_rev", ""),
                    "evalue":     row.get("evalue", ""),
                }
    return rbh


def find_orthofinder_pairwise(of_dir: str, species: str = "Aniger_iJB1325") -> str | None:
    """Locate the <species> vs Fvenenatum pairwise orthologs file."""
    pattern = os.path.join(
        of_dir, "**", "Orthologues",
        f"Orthologues_{species}",
        f"{species}__v__Fvenenatum.tsv",
    )
    matches = glob.glob(pattern, recursive=True)
    return matches[0] if matches else None


def load_orthofinder(path: str) -> dict[str, list[str]]:
    """Parse OrthoFinder pairwise orthologs TSV.

    Format: Orthogroup TAB Aniger_genes(comma) TAB Fvenenatum_genes(comma)
    Returns {iJB1325_gene: [F.venenatum_genes...]}.
    """
    out: dict[str, list[str]] = defaultdict(list)
    with open(path) as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            _og, aniger_field, fven_field = parts[0], parts[1], parts[2]
            aniger_genes = [g.strip() for g in aniger_field.split(",") if g.strip()]
            fven_genes = [g.strip() for g in fven_field.split(",") if g.strip()]
            for ag in aniger_genes:
                out[ag].extend(fven_genes)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--rbh-tsv",     default="orthology/orthologs_cbs513_rbh.tsv",
                        help="RBH ortholog table from build_orthologs.py (will be backed up)")
    parser.add_argument("--of-dir",      default="orthology/orthofinder_output",
                        help="OrthoFinder output root directory")
    parser.add_argument("--output",      default="orthology/orthologs_cbs513.tsv",
                        help="Merged output")
    parser.add_argument("--aniger-species", default="Aniger_CBS513",
                        help="OrthoFinder species name of the A. niger input file")
    parser.add_argument("--backup",      default="orthology/orthologs_cbs513_rbh_only.tsv",
                        help="Where to save the original RBH-only table")
    args = parser.parse_args()

    if not os.path.isfile(args.rbh_tsv):
        sys.exit(f"ERROR: RBH ortholog table not found: {args.rbh_tsv}")
    of_pairwise = find_orthofinder_pairwise(args.of_dir, args.aniger_species)
    if not of_pairwise:
        sys.exit(f"ERROR: OrthoFinder pairwise file not found under {args.of_dir}")
    print(f"OrthoFinder pairwise file: {of_pairwise}")

    # Backup the RBH-only table once
    if not os.path.isfile(args.backup):
        os.system(f"cp {args.rbh_tsv} {args.backup}")
        print(f"Backed up RBH-only table to {args.backup}")

    rbh = load_rbh(args.rbh_tsv)
    of  = load_orthofinder(of_pairwise)
    print(f"RBH orthologs loaded:        {len(rbh)}")
    print(f"OrthoFinder iJB1325 genes:   {len(of)}")

    rbh_keys = set(rbh.keys())
    of_keys  = set(of.keys())
    both     = rbh_keys & of_keys
    rbh_only = rbh_keys - of_keys
    of_only  = of_keys  - rbh_keys

    # Sanity check: where both agree, do they map to the same F. venenatum gene?
    agree = 0
    disagree = 0
    for g in both:
        if rbh[g]["target"] in of[g]:
            agree += 1
        else:
            disagree += 1

    print()
    print(f"Source breakdown:")
    print(f"  RBH only:                  {len(rbh_only):>5}  (~lost in OF)")
    print(f"  OrthoFinder only:          {len(of_only):>5}  (NEW augmenting RBH)")
    print(f"  Both (intersection):       {len(both):>5}")
    print(f"    same target:             {agree:>5}  (high-confidence)")
    print(f"    different target:        {disagree:>5}  (RBH wins; flag for review)")

    # Merge: for each template gene, write one row
    rows: list[dict[str, str]] = []
    all_template = rbh_keys | of_keys
    for tg in sorted(all_template):
        if tg in rbh:
            target = rbh[tg]["target"]
            source = "RBH+OF" if tg in of and target in of[tg] else (
                "RBH"   if tg not in of else "RBH(OF-disagree)")
            row = {
                "template_gene": tg,
                "target_gene":   target,
                "pident_fwd":    rbh[tg].get("pident_fwd", ""),
                "pident_rev":    rbh[tg].get("pident_rev", ""),
                "evalue":        rbh[tg].get("evalue", ""),
                "method":        source,
            }
        else:
            # OrthoFinder-only: pick the first F. venenatum candidate
            target = of[tg][0]
            row = {
                "template_gene": tg,
                "target_gene":   target,
                "pident_fwd":    "",
                "pident_rev":    "",
                "evalue":        "",
                "method":        "OF",
            }
        rows.append(row)

    fields = ["template_gene", "target_gene", "pident_fwd", "pident_rev", "evalue", "method"]
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print()
    print(f"Merged orthologs written: {args.output}")
    print(f"  total template genes with ortholog: {len(rows)}")
    src_counter = defaultdict(int)
    for r in rows:
        src_counter[r["method"]] += 1
    for src, n in sorted(src_counter.items()):
        print(f"    {src:20s} {n:>5}")


if __name__ == "__main__":
    main()
