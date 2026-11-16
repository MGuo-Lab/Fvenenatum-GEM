"""Build the released model iJL1143 from the Curated model.

Only the gene layer changes. Reactions, metabolites, bounds and the objective stay as they are, and the
script stops if any of them differ.

The Curated SBML holds 1,203 gene objects. 1,142 are F. venenatum genes. The other 61 are:
  - 52 template genes that no GPR references. They are the genes of the template reactions with no ortholog.
    curate_model.needs_based_gapfill adds those reactions temporarily, and clearing their GPRs leaves the
    gene objects behind.
  - 9 placeholders, OLD_<template gene>, for complex subunits with no ortholog in the orthology step.
    RAVEN uses this form in draft models for genes "not belonging to the new organism".

Each placeholder was aligned to the whole proteome (check_unmapped_subunits.py). The carrier subunit of
molybdopterin synthase has a homologue of the same function, so it is replaced. The other eight are dropped
from their GPRs. Every change is written to model/curation/release_gene_changes.csv.

Usage:  python build_release_model.py
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import cobra
from cobra.manipulation import remove_genes

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "build" / "curated.xml"                 # written by curate_model.py
PROTEOME = HERE / "raw_data" / "fven_proteins.faa"
OUTPUT = HERE.parent / "model" / "iJL1143.xml"
CHANGES = HERE.parent / "model" / "curation" / "release_gene_changes.csv"
MODEL_ID = "iJL1143"
MODEL_NAME = "iJL1143: genome-scale metabolic model of Fusarium venenatum A3/5, version 1.1.0"

# The nine placeholders: replacement gene (None = dropped) and the reason.
SUBUNITS = {
    "OLD_An12g01120": ("A0A2L2TNV3", "molybdopterin synthase sulfur carrier subunit; the hit has the same UniProt "
                                     "function and meets the thresholds, and it pairs with paralogue An13g01070, "
                                     "which is not in the template"),
    "OLD_An16g04960": (None, "fragment of the carbamoyl phosphate synthase large chain; the hit covers all of it "
                             "and is already in the GPR through An16g04970"),
    "OLD_An09g03020": (None, "the template rule is (subunit and An03g02930) or An03g02930, so the rule reduces to "
                             "A0A2L2T4R4 whatever the subunit maps to"),
    "OLD_An05g02310": (None, "below the thresholds; nearest protein is SUR1, already in the GPR"),
    "OLD_An09g00100": (None, "below the thresholds; nearest protein is SUR1, already in the GPR"),
    "OLD_An09g01820": (None, "A. niger secondary metabolism; below the thresholds, family homologue only"),
    "OLD_An09g01860": (None, "A. niger secondary metabolism; below the thresholds, family homologue only"),
    "OLD_An10g00180": (None, "A. niger secondary metabolism; below the thresholds"),
    "OLD_An11g07330": (None, "A. niger secondary metabolism; below the thresholds, family homologue only"),
}

# The resulting rules, written out so that each can be read and checked.
RULES = {
    "r479": "A0A2L2T2T6 and A0A2L2SWY4",
    "r856": "A0A2L2TRN1",
    "r857": "A0A2L2TRN1",
    "r1953": "(A0A2L2T9N6 or A0A2L2TDU4) and A0A2L2T7V9",
    "r1997": "A0A2L2U5I3",
    "r1998": "(A0A2L2SWF0 or A0A2L2SYG4 or A0A2L2TCR7 or A0A2L2TGP7 or A0A2L2TK98 or A0A2L2TVL6) and A0A2L2TDQ1",
    "r2003": "A0A2L2TQR4",
    "r2328": "A0A2L2TK06 and A0A2L2T0A5 and A0A2L2TNV3",
    "r2329": "A0A2L2T4R4",
}


def network(model: cobra.Model) -> dict:
    """Everything that defines the flux problem. Gene associations are left out on purpose."""
    return {
        "reactions": [(r.id, r.name, r.lower_bound, r.upper_bound, r.objective_coefficient,
                       sorted((m.id, c) for m, c in r.metabolites.items())) for r in model.reactions],
        "metabolites": [(m.id, m.name, m.formula, m.charge, m.compartment) for m in model.metabolites],
    }


def main() -> None:
    model = cobra.io.read_sbml_model(str(SOURCE))
    before = network(model)
    proteome = set(re.findall(r"^>\w\w\|([^|]+)\|", PROTEOME.read_text(encoding="utf-8"), re.M))
    old_rules = {r.id: r.gene_reaction_rule for r in model.reactions}
    rows = []

    for rid, rule in RULES.items():
        reaction = model.reactions.get_by_id(rid)
        placeholders = {g.id for g in reaction.genes if g.id.startswith("OLD_")}
        kept = {g.id for g in reaction.genes} - placeholders
        reaction.gene_reaction_rule = rule
        added = {SUBUNITS[p][0] for p in placeholders if SUBUNITS[p][0]}
        if {g.id for g in reaction.genes} != kept | added:
            raise ValueError(f"{rid}: new rule does not match the old genes plus replacements")

    for gene_id, (replacement, reason) in SUBUNITS.items():
        rids = sorted(r for r in RULES if gene_id in old_rules[r])
        rows.append({"gene": gene_id, "kind": "unmapped complex subunit",
                     "action": f"replaced by {replacement}" if replacement else "dropped from GPR",
                     "reactions": " ".join(rids), "reason": reason})

    orphans = [g for g in model.genes if not g.reactions]
    for gene in (g for g in orphans if not g.id.startswith("OLD_")):
        rows.append({"gene": gene.id, "kind": "template gene in no GPR", "action": "removed",
                     "reactions": "", "reason": "gene of a template reaction with no ortholog"})
    remove_genes(model, orphans, remove_reactions=False)

    new_gene = model.genes.get_by_id("A0A2L2TNV3")
    new_gene.name = "G_A0A2L2TNV3"

    if network(model) != before:
        raise ValueError("the metabolic network changed")
    foreign = [g.id for g in model.genes if g.id not in proteome]
    if foreign:
        raise ValueError(f"genes outside the proteome: {foreign}")
    if any(not g.reactions for g in model.genes):
        raise ValueError("unreferenced genes remain")
    changed = sorted(r.id for r in model.reactions if r.gene_reaction_rule != old_rules[r.id])
    if changed != sorted(RULES):
        raise ValueError(f"unexpected GPR changes: {changed}")

    model.id, model.name = MODEL_ID, MODEL_NAME
    cobra.io.write_sbml_model(model, str(OUTPUT))
    with open(CHANGES, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"{MODEL_ID}: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites, "
          f"{len(model.genes)} genes; {len(rows)} gene changes written to {CHANGES.name}")
    if len(model.genes) != int(MODEL_ID[3:]):
        raise ValueError("the gene count no longer matches the model name")


if __name__ == "__main__":
    main()
