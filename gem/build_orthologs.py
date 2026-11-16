"""Map template genes to F. venenatum by reciprocal best hits (RBH) from blastp.

Template: the CBS 513.88 version of iJB1325 (Brandl et al. 2018, Additional file 3), whose gene ids are A. niger
CBS 513.88 locus tags (An01g00010, ...). They resolve to the UniProt proteome UP000006706.

  1. Write one A. niger sequence per template gene, with the template gene id as FASTA id.
  2. blastp the template sequences against the F. venenatum A3/5 proteome (UniProt UP000245910), and the whole
     F. venenatum proteome against the whole A. niger proteome.
  3. Keep the best hit per query that passes the identity and query-coverage thresholds (40 %, 70 %,
     E <= 1e-30), and keep the pairs that are each other's best hit.
  4. Write template_gene, target_gene, pident_fwd, pident_rev, evalue.

blastp and makeblastdb (BLAST+) must be on PATH. The orthogroups that complement these pairs come from OrthoFinder
(jobs/orthology.sbatch); merge_orthologs.py combines the two.

Usage:
    python build_orthologs.py --output orthology/orthologs_rbh.tsv --threads 16
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys

import cobra


HERE = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(HERE, "raw_data")
ORTHO_DIR = os.path.join(HERE, "orthology")
BLAST_DIR = os.path.join(ORTHO_DIR, "blast")

DEFAULT_TEMPLATE = os.path.join(RAW_DIR, "MOESM3_cbs51388.xml")
DEFAULT_ANIGER_FA = os.path.join(RAW_DIR, "aniger_cbs513.faa")
DEFAULT_TARGET_FA = os.path.join(RAW_DIR, "fven_proteins.faa")
DEFAULT_OUTPUT = os.path.join(ORTHO_DIR, "orthologs_rbh.tsv")


# ── Template -> A. niger representative proteome ────────────────────────────

_HEADER_RE = re.compile(r">(\S+)")


def parse_aniger_fasta(path: str, source: str = "direct") -> dict[str, tuple[str, str]]:
    """Return {gene_id: (header_id, full_sequence)}, keyed by the first token of each FASTA header."""
    entries: dict[str, tuple[str, str]] = {}
    cur_id = None
    cur_header = None
    cur_seq: list[str] = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur_id is not None:
                    entries[cur_id] = (cur_header, "".join(cur_seq))
                m_h = _HEADER_RE.match(line)
                cur_id = cur_header = m_h.group(1) if m_h else None
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_id is not None:
            entries[cur_id] = (cur_header, "".join(cur_seq))
    return entries


def build_template_proteome(
    template_model_path: str,
    aniger_fa_path: str,
    aniger_source: str,
    out_fa_path: str,
) -> tuple[int, int]:
    """Write a FASTA containing one A. niger sequence per recoverable iJB1325 gene.

    Headers are formatted as ``>{template_gene_id}|{header_id}`` so the BLAST
    output can be parsed back to template gene IDs without a separate mapping
    table.

    Returns (template_gene_count, recovered_count).
    """
    print(f"Loading template model: {template_model_path}")
    model = cobra.io.read_sbml_model(template_model_path)
    template_gene_count = len(model.genes)
    print(f"  {template_gene_count} template genes")

    print(f"Parsing A. niger proteome: {aniger_fa_path} (source={aniger_source})")
    aspn_db = parse_aniger_fasta(aniger_fa_path, aniger_source)
    print(f"  {len(aspn_db):,} entries indexed by ID")

    recovered = 0
    src_via_id = 0
    src_via_name = 0
    os.makedirs(os.path.dirname(out_fa_path), exist_ok=True)
    with open(out_fa_path, "w") as out:
        for g in model.genes:
            gid = g.id
            gname = (g.name or "").strip()
            seq_entry = None
            via = None
            if gid in aspn_db:
                seq_entry = aspn_db[gid]
                via = "id"
            elif gname.isdigit() and gname in aspn_db:
                seq_entry = aspn_db[gname]
                via = "name"
            if seq_entry is None:
                continue
            header_id, seq = seq_entry
            # the template gene id alone is the FASTA tag (BLAST limits local ids to 50 characters)
            tag = gid
            out.write(f">{tag}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
            recovered += 1
            if via == "id":
                src_via_id += 1
            else:
                src_via_name += 1

    print(f"\nTemplate proteome -> {out_fa_path}")
    print(f"  recovered template genes: {recovered}/{template_gene_count} "
          f"({recovered / template_gene_count * 100:.1f}%)")
    print(f"  via gene.id: {src_via_id}")
    if src_via_name:
        print(f"  via gene.name: {src_via_name}")
    return template_gene_count, recovered


# ── BLAST RBH ───────────────────────────────────────────────────────────────

BLAST_FIELDS = (
    "qseqid sseqid pident length mismatch gapopen qstart qend "
    "sstart send evalue bitscore qlen slen"
)


def run_blast(
    query: str,
    db_prefix: str,
    out_path: str,
    evalue: float,
    threads: int,
) -> None:
    """Run blastp with the configured tabular output format."""
    print(f"  blastp {os.path.basename(query)} -> "
          f"{os.path.basename(db_prefix)} (E<={evalue}, threads={threads})")
    cmd = [
        "blastp",
        "-query", query,
        "-db", db_prefix,
        "-out", out_path,
        "-evalue", str(evalue),
        "-outfmt", f"6 {BLAST_FIELDS}",
        "-num_threads", str(threads),
        "-max_target_seqs", "5",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"ERROR: blastp failed:\n{res.stderr}")
    if res.stderr.strip():
        # BLAST often emits informational warnings on stderr; surface them.
        print("    [blast stderr]:", res.stderr.strip().splitlines()[0][:200])


def make_blast_db(fasta: str, db_prefix: str) -> None:
    print(f"  makeblastdb {os.path.basename(fasta)} -> {os.path.basename(db_prefix)}")
    cmd = [
        "makeblastdb",
        "-in", fasta,
        "-dbtype", "prot",
        "-out", db_prefix,
        "-parse_seqids",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"ERROR: makeblastdb failed:\n{res.stderr}")


def best_hit_per_query(
    blast_out: str,
    pident_min: float,
    coverage_min: float,
) -> dict[str, tuple[str, float, float, float]]:
    """Return {query_id: (subject_id, pident, evalue, bitscore)} keeping the
    highest-bitscore hit that passes identity & query-coverage thresholds.
    """
    best: dict[str, tuple[str, float, float, float]] = {}
    with open(blast_out) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 14:
                continue
            qid, sid = parts[0], parts[1]
            pident = float(parts[2])
            length = int(parts[3])
            evalue = float(parts[10])
            bitscore = float(parts[11])
            qlen = int(parts[12])
            qcov = (length / qlen) * 100.0 if qlen > 0 else 0.0
            if pident < pident_min or qcov < coverage_min:
                continue
            cur = best.get(qid)
            if cur is None or bitscore > cur[3]:
                best[qid] = (sid, pident, evalue, bitscore)
    return best


def parse_template_gene_from_query(qid: str) -> str:
    """Recover the template gene ID from FASTA headers we wrote earlier.

    Header form: ``{template_gene}``, which keeps the BLAST local id within the 50-character limit.
    """
    return qid.split("|", 1)[0]


_UNIPROT_HEADER_RE = re.compile(r"^(?:sp|tr)\|([A-Z0-9]+)\|")


def clean_target_id(sseqid: str) -> str:
    """Extract a SBML-safe identifier from a target proteome BLAST hit.

    The F. venenatum FASTA uses UniProt headers like ``tr|A0A2L2SSM9|A0A2L2SSM9_9HYPO``.
    SBML/COBRApy GPR parsing requires gene IDs to be valid Python identifiers, so
    pipe characters and the db/entry-name fields are stripped down to the
    UniProt accession (``A0A2L2SSM9``). For non-UniProt headers, fall back to
    replacing illegal characters with underscores.
    """
    m = _UNIPROT_HEADER_RE.match(sseqid)
    if m:
        return m.group(1)
    return re.sub(r"[^A-Za-z0-9_]", "_", sseqid)


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--template-model", default=DEFAULT_TEMPLATE)
    parser.add_argument("--target-fasta",   default=DEFAULT_TARGET_FA,
                        help="F. venenatum proteome FASTA")
    parser.add_argument("--output",         default=DEFAULT_OUTPUT,
                        help="orthologs.tsv path")
    parser.add_argument("--evalue",         type=float, default=1e-30,
                        help="Maximum BLAST e-value (default: 1e-30, "
                             "stricter than the canonical 1e-10 to reduce "
                             "twilight-zone false-positive orthologs).")
    parser.add_argument("--pident",         type=float, default=40.0,
                        help="Minimum percent identity (default: 40, the "
                             "twilight-zone boundary — orthologs below 30%% "
                             "are unreliable for closely related fungi).")
    parser.add_argument("--coverage",       type=float, default=70.0,
                        help="Minimum query coverage %% (default: 70, in line "
                             "with publication-grade RBH practice).")
    parser.add_argument("--threads",        type=int, default=os.cpu_count() or 4)
    parser.add_argument("--keep-blast",     action="store_true",
                        help="Keep intermediate BLAST tables and DBs")
    parser.add_argument("--aniger-fasta",   default=DEFAULT_ANIGER_FA,
                        help="A. niger proteome whose FASTA ids are the template "
                             "gene ids (CBS 513.88 An locus tags). The reverse "
                             "search runs against this whole proteome, as "
                             "standard RBH does.")
    args = parser.parse_args()

    if shutil.which("blastp") is None or shutil.which("makeblastdb") is None:
        sys.exit("ERROR: blastp/makeblastdb (BLAST+) not on PATH.")

    if not os.path.isfile(args.target_fasta):
        sys.exit(f"ERROR: target proteome not found: {args.target_fasta}")
    if not os.path.isfile(args.aniger_fasta):
        sys.exit(f"ERROR: A. niger proteome not found: {args.aniger_fasta}")
    if not os.path.isfile(args.template_model):
        sys.exit(f"ERROR: template model not found: {args.template_model}")

    os.makedirs(BLAST_DIR, exist_ok=True)

    # 1. A. niger proteome, FASTA ids = template gene ids
    aniger_fa, aniger_source = args.aniger_fasta, "direct"

    # 2. Build template-representative proteome
    template_fa = os.path.join(BLAST_DIR, "template_iJB1325_repr.faa")
    n_template, n_recovered = build_template_proteome(
        args.template_model, aniger_fa, aniger_source, template_fa,
    )
    if n_recovered == 0:
        sys.exit("ERROR: no template genes could be matched to the A. niger proteome.")

    # 3. BLAST databases
    print("\nBuilding BLAST databases...")
    db_template = os.path.join(BLAST_DIR, "db_template")
    db_target   = os.path.join(BLAST_DIR, "db_target")
    # the reverse search sees every A. niger protein, so a target whose true best hit is a non-template paralog is
    # not forced onto a template gene
    make_blast_db(aniger_fa, db_template)
    make_blast_db(args.target_fasta,  db_target)

    # 4. Reciprocal blastp
    print("\nRunning blastp (forward: template -> target)...")
    fwd = os.path.join(BLAST_DIR, "fwd.tsv")
    run_blast(template_fa, db_target, fwd, args.evalue, args.threads)

    print("Running blastp (reverse: target -> template)...")
    rev = os.path.join(BLAST_DIR, "rev.tsv")
    run_blast(args.target_fasta, db_template, rev, args.evalue, args.threads)

    # 5. Best-hit per query, then RBH intersection
    print("\nIdentifying reciprocal best hits...")
    fwd_best = best_hit_per_query(fwd, args.pident, args.coverage)
    rev_best = best_hit_per_query(rev, args.pident, args.coverage)
    print(f"  forward queries with hit: {len(fwd_best):,}")
    print(f"  reverse queries with hit: {len(rev_best):,}")

    rbh: list[tuple[str, str, float, float, float]] = []
    for q_temp, (s_target, pid, ev, bs) in fwd_best.items():
        rev_hit = rev_best.get(s_target)
        if rev_hit is None:
            continue
        if rev_hit[0] != q_temp:
            continue
        template_gene = parse_template_gene_from_query(q_temp)
        target_gene = clean_target_id(s_target)
        # Average pident from both directions for the report
        rbh.append((template_gene, target_gene, pid, rev_hit[1], ev))

    # Multiple template genes can share an A. niger protein (we tagged each).
    # Each row in our forward results corresponds to one template gene.
    print(f"  reciprocal best hit pairs: {len(rbh):,}")

    # 6. Write orthologs.tsv
    seen: set[str] = set()
    rows: list[dict[str, str]] = []
    for tpl, tgt, pid_fwd, pid_rev, evalue in rbh:
        if tpl in seen:
            continue
        seen.add(tpl)
        rows.append({
            "template_gene": tpl,
            "target_gene": tgt,
            "pident_fwd": f"{pid_fwd:.2f}",
            "pident_rev": f"{pid_rev:.2f}",
            "evalue":      f"{evalue:.2e}",
        })

    rows.sort(key=lambda r: r["template_gene"])
    fields = ["template_gene", "target_gene", "pident_fwd", "pident_rev", "evalue"]
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"\nOrthologs written: {args.output}")
    print(f"  template genes total:           {n_template}")
    print(f"  template genes with sequence:   {n_recovered}  "
          f"({n_recovered / n_template * 100:.1f}%)")
    print(f"  template genes with RBH:        {len(rows)}  "
          f"({len(rows) / n_template * 100:.1f}%)")

    if not args.keep_blast:
        shutil.rmtree(BLAST_DIR, ignore_errors=True)
        print(f"  cleaned intermediate BLAST files in {BLAST_DIR}")


if __name__ == "__main__":
    main()
