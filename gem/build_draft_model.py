"""Build a draft genome-scale model from a template and ortholog mapping.

Steps:
  1. Load a template SBML model.
  2. Load an ortholog mapping table (TSV/CSV) that maps template gene IDs
     to target organism gene IDs.
  3. Replace gene identifiers in GPR rules using the mapping.
  4. Remove reactions whose GPR genes could not be mapped (no orthologs).
  5. Remove orphan metabolites and genes.
  6. Save the draft model (build/draft.xml by default).

Ortholog mapping file format (TSV, with header):
    template_gene<TAB>target_gene
    b0001<TAB>A35_00001
    b0002<TAB>A35_00002
    ...

Rows where target_gene is empty are treated as unmapped.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from copy import deepcopy

import cobra


DEFAULT_OUTPUT = os.path.join(os.path.dirname(__file__), "build", "draft.xml")


# ── Ortholog mapping ─────────────────────────────────────────────────────────

def load_ortholog_mapping(path: str) -> dict[str, str]:
    """Load ortholog mapping from a TSV/CSV file.

    Returns a dict mapping template_gene -> target_gene.
    """
    mapping: dict[str, str] = {}

    with open(path, newline="") as f:
        # Auto-detect delimiter
        sample = f.read(4096)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        reader = csv.DictReader(f, dialect=dialect)

        # Identify columns: expect 'template_gene' and 'target_gene'
        # or fall back to first two columns
        fieldnames = reader.fieldnames
        if fieldnames is None or len(fieldnames) < 2:
            sys.exit(f"ERROR: Mapping file must have at least 2 columns: {path}")

        col_template = None
        col_target = None
        for name in fieldnames:
            lower = name.strip().lower().replace(" ", "_")
            if "template" in lower or lower in ("source", "from", "query"):
                col_template = name
            elif "target" in lower or lower in ("subject", "to", "hit"):
                col_target = name

        if col_template is None or col_target is None:
            # Fall back to first two columns
            col_template = fieldnames[0]
            col_target = fieldnames[1]

        for row in reader:
            tpl = row[col_template].strip()
            tgt = row[col_target].strip()
            if tpl and tgt:
                mapping[tpl] = tgt

    return mapping


# ── GPR rewriting ────────────────────────────────────────────────────────────
#
# GPRs are parsed into a small Boolean tree and remapped node by node, so the
# output is always well-formed (the previous regex/strip approach could emit
# unbalanced parentheses, which cobra silently turns into an empty GPR).

_Node = tuple  # ("gene", id) | ("and", [nodes]) | ("or", [nodes])


def _parse_gpr(gpr_string: str) -> _Node:
    """Parse 'a and (b or c)' into a Boolean tree (and binds tighter than or)."""
    tokens = re.findall(r"\(|\)|[^\s()]+", gpr_string)
    pos = 0

    def peek() -> str | None:
        return tokens[pos] if pos < len(tokens) else None

    def take() -> str:
        nonlocal pos
        pos += 1
        return tokens[pos - 1]

    def expr() -> _Node:
        terms = [term()]
        while (peek() or "").lower() == "or":
            take()
            terms.append(term())
        return terms[0] if len(terms) == 1 else ("or", terms)

    def term() -> _Node:
        factors = [factor()]
        while (peek() or "").lower() == "and":
            take()
            factors.append(factor())
        return factors[0] if len(factors) == 1 else ("and", factors)

    def factor() -> _Node:
        tok = take()
        if tok == "(":
            node = expr()
            if take() != ")":
                raise ValueError(f"Unbalanced GPR: {gpr_string!r}")
            return node
        if tok in (")",) or tok.lower() in ("and", "or"):
            raise ValueError(f"Malformed GPR: {gpr_string!r}")
        return ("gene", tok)

    tree = expr()
    if pos != len(tokens):
        raise ValueError(f"Trailing tokens in GPR: {gpr_string!r}")
    return tree


def _genes(node: _Node) -> list[str]:
    if node[0] == "gene":
        return [node[1]]
    return [g for child in node[1] for g in _genes(child)]


def _to_string(node: _Node, parent: str | None = None) -> str:
    if node[0] == "gene":
        return node[1]
    op = node[0]
    s = f" {op} ".join(_to_string(c, op) for c in node[1])
    return f"({s})" if parent is not None and parent != op else s


def remap_gpr_detailed(
    gpr_string: str,
    mapping: dict[str, str],
    complex_policy: str = "strict",
    fallback: dict[str, list[str]] | None = None,
) -> tuple[str | None, dict]:
    """Rewrite a GPR from template genes to target genes.

    Evidence tiers per template gene:
      tier1 — 1:1 ortholog from ``mapping``;
      tier2 — F. venenatum members of the same orthogroup (``fallback``),
              used only when there is no 1:1 ortholog.

    'or' (isozymes): unmapped members are dropped; None if none remain.
    'and' (complex):
      complex_policy="strict" — None if any subunit is unmapped;
      complex_policy="raven"  — kept if at least one subunit maps, unmapped
        subunits recorded as ``OLD_<template gene>`` (RAVEN
        getModelFromHomology convention).

    Returns (new_gpr_or_None, info) with info = {"tiers": set,
    "complex_partial": bool}. An empty rule returns "".
    """
    if complex_policy not in ("strict", "raven"):
        raise ValueError(f"Unknown complex_policy: {complex_policy}")
    info = {"tiers": set(), "complex_partial": False}
    if not gpr_string or not gpr_string.strip():
        return "", info
    fallback = fallback or {}

    def walk(node: _Node) -> _Node | None:
        if node[0] == "gene":
            gene = node[1]
            if gene in mapping:
                info["tiers"].add("tier1")
                return ("gene", mapping[gene])
            if fallback.get(gene):
                info["tiers"].add("tier2")
                targets = sorted(fallback[gene])
                if len(targets) == 1:
                    return ("gene", targets[0])
                return ("or", [("gene", t) for t in targets])
            return None
        children = [walk(c) for c in node[1]]
        if node[0] == "or":
            kept = [c for c in children if c is not None]
            if not kept:
                return None
            return kept[0] if len(kept) == 1 else ("or", kept)
        # and-node
        if all(c is not None for c in children):
            return ("and", children)
        if complex_policy == "raven" and any(c is not None for c in children):
            info["complex_partial"] = True
            filled = []
            for orig, new in zip(node[1], children):
                if new is not None:
                    filled.append(new)
                else:
                    olds = [("gene", f"OLD_{g}") for g in _genes(orig)]
                    filled.append(olds[0] if len(olds) == 1 else ("or", olds))
            return ("and", filled)
        return None

    remapped = walk(_parse_gpr(gpr_string))
    return (None if remapped is None else _to_string(remapped)), info


def remap_gpr(
    gpr_string: str,
    mapping: dict[str, str],
    complex_policy: str = "strict",
    fallback: dict[str, list[str]] | None = None,
) -> str | None:
    """Rewrite a GPR rule; see :func:`remap_gpr_detailed`."""
    return remap_gpr_detailed(gpr_string, mapping, complex_policy, fallback)[0]


def load_orthogroups(path: str, template_col: str = "Aniger_iJB1325",
                     target_col: str = "Fvenenatum") -> dict[str, list[str]]:
    """Read OrthoFinder ``Orthogroups.tsv`` into {template_gene: [target genes]}.

    Only template genes whose orthogroup contains at least one target gene
    are returned.
    """
    out: dict[str, set[str]] = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            tpl = [g.strip() for g in (row.get(template_col) or "").split(",") if g.strip()]
            tgt = [g.strip() for g in (row.get(target_col) or "").split(",") if g.strip()]
            if not tgt:
                continue
            for g in tpl:
                out.setdefault(g, set()).update(tgt)
    return {g: sorted(v) for g, v in out.items()}


# ── Model building ───────────────────────────────────────────────────────────

def build_draft_model(
    template: cobra.Model,
    mapping: dict[str, str],
    keep_unmapped: bool = False,
    complex_policy: str = "strict",
    fallback: dict[str, list[str]] | None = None,
    evidence_out: list[dict] | None = None,
) -> cobra.Model:
    """Build a draft model by remapping GPRs.

    ``complex_policy`` and ``fallback`` are passed to :func:`remap_gpr_detailed`.
    If ``evidence_out`` is a list, one row per template reaction is appended:
    {reaction_id, template_gpr, new_gpr, status} with status one of
    tier1 / tier2_orthogroup / complex_partial / no_gpr_in_template /
    stripped (lenient) / removed (strict).

    With ``keep_unmapped=False`` (strict mode, the default) reactions whose GPR
    cannot be mapped to a target ortholog are removed. This is biologically
    rigorous but requires high orthology coverage; with low coverage the
    resulting draft will be missing many essential pathways and will not grow
    on minimal medium.

    With ``keep_unmapped=True`` (lenient/homology mode) the same reactions are
    kept but their GPR is stripped. The metabolic network topology is fully
    preserved, at the cost of losing gene-level evidence for those reactions.
    Gene-knockout analyses on those reactions become uninformative.
    """

    draft = deepcopy(template)
    # Template's id is empty; set a meaningful identifier on the draft
    draft.id = "fven_a3_5_draft"
    draft.name = "Fusarium venenatum A3/5 draft GEM (from iJB1325)"

    reactions_to_remove = []
    mapped_count = 0
    kept_no_gpr = 0
    stripped_count = 0

    for rxn in draft.reactions:
        gpr_str = rxn.gene_reaction_rule.strip()

        # Keep exchange, demand, sink, and spontaneous reactions (no GPR)
        if not gpr_str:
            kept_no_gpr += 1
            status, new_gpr = "no_gpr_in_template", ""
        else:
            new_gpr, info = remap_gpr_detailed(gpr_str, mapping, complex_policy, fallback)
            if new_gpr is None:
                if keep_unmapped:
                    # Lenient mode: keep the reaction, strip the GPR.
                    rxn.gene_reaction_rule = ""
                    stripped_count += 1
                    status = "stripped"
                else:
                    reactions_to_remove.append(rxn)
                    status = "removed"
            else:
                rxn.gene_reaction_rule = new_gpr
                mapped_count += 1
                if info["complex_partial"]:
                    status = "complex_partial"
                elif "tier2" in info["tiers"]:
                    status = "tier2_orthogroup"
                else:
                    status = "tier1"
        if evidence_out is not None:
            evidence_out.append({"reaction_id": rxn.id, "template_gpr": gpr_str,
                                 "new_gpr": new_gpr or "", "status": status})

    if reactions_to_remove:
        draft.remove_reactions(reactions_to_remove, remove_orphans=True)

    # Prune orphan genes inherited from the template's <fbc:listOfGeneProducts>.
    # In lenient mode, stripping GPRs leaves the gene definitions in SBML even
    # though no reaction references them. This inflates the gene count without
    # affecting FBA — clean it up so model.genes reflects reality.
    orphan_genes = [g for g in draft.genes if not g.reactions]
    if orphan_genes:
        cobra.manipulation.delete.remove_genes(draft, orphan_genes, remove_reactions=False)

    # Summary
    print(f"\nDraft model summary:")
    print(f"  Template reactions:           {len(template.reactions)}")
    print(f"  Reactions with GPR mapped:    {mapped_count}")
    print(f"  Reactions without GPR kept:   {kept_no_gpr}")
    if keep_unmapped:
        print(f"  Reactions kept w/ stripped GPR: {stripped_count}  (lenient mode)")
    else:
        print(f"  Reactions removed (unmapped): {len(reactions_to_remove)}")
    print(f"  Draft reactions:              {len(draft.reactions)}")
    print(f"  Draft metabolites:            {len(draft.metabolites)}")
    print(f"  Draft genes (linked to rxns): {len(draft.genes)}")
    if orphan_genes:
        print(f"  Orphan template-genes pruned: {len(orphan_genes)}")

    return draft


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a draft metabolic model from a template + ortholog mapping"
    )
    parser.add_argument("template", help="Path to template SBML model")
    parser.add_argument("mapping", help="Path to ortholog mapping file (TSV/CSV)")
    parser.add_argument("--output", default=DEFAULT_OUTPUT,
                        help="Output SBML path (default: build/draft.xml)")
    parser.add_argument("--keep-unmapped", action="store_true",
                        help="Lenient/homology mode: keep reactions whose GPR cannot "
                             "be remapped, just strip the GPR. Trades gene-level "
                             "evidence for full metabolic network coverage. Use this "
                             "when ortholog coverage is too low for strict mode to "
                             "produce a functional model.")
    parser.add_argument("--complex-policy", choices=("strict", "raven"), default="strict",
                        help="strict: drop an enzyme complex if any subunit is unmapped; "
                             "raven: keep it if at least one subunit maps, recording "
                             "unmapped subunits as OLD_<gene> (RAVEN convention)")
    parser.add_argument("--orthogroups", default=None,
                        help="OrthoFinder Orthogroups.tsv; template genes without a 1:1 "
                             "ortholog fall back to same-orthogroup target genes (tier 2)")
    parser.add_argument("--evidence", default=None,
                        help="Write the per-reaction evidence table (TSV) here")
    args = parser.parse_args()

    # Load template
    if not os.path.isfile(args.template):
        sys.exit(f"ERROR: Template model not found: {args.template}")
    print(f"Loading template: {args.template}")
    template = cobra.io.read_sbml_model(args.template)
    print(f"  {template.id}: {len(template.reactions)} reactions, "
          f"{len(template.metabolites)} metabolites, {len(template.genes)} genes")

    # Load ortholog mapping
    if not os.path.isfile(args.mapping):
        sys.exit(f"ERROR: Mapping file not found: {args.mapping}")
    print(f"Loading ortholog mapping: {args.mapping}")
    mapping = load_ortholog_mapping(args.mapping)
    print(f"  {len(mapping)} gene mappings loaded")

    fallback = None
    if args.orthogroups:
        fallback = load_orthogroups(args.orthogroups)
        print(f"  {len(fallback)} template genes with orthogroup fallback")

    # Build draft
    evidence: list[dict] = []
    draft = build_draft_model(template, mapping, keep_unmapped=args.keep_unmapped,
                              complex_policy=args.complex_policy, fallback=fallback,
                              evidence_out=evidence)
    if args.evidence:
        os.makedirs(os.path.dirname(os.path.abspath(args.evidence)), exist_ok=True)
        with open(args.evidence, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(evidence[0]), delimiter="\t")
            w.writeheader()
            w.writerows(evidence)
        print(f"Evidence table saved to: {args.evidence}")

    # Save
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    cobra.io.write_sbml_model(draft, args.output)
    print(f"\nDraft model saved to: {args.output}")


if __name__ == "__main__":
    main()
