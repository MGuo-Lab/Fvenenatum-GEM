"""Data-independent curation of the F. venenatum GEM (strict-homology path).

Builds the Curated model from the iJB1325 template (CBS 513.88 version) without
looking at any fermentation data. build_release_model.py then cleans its gene
layer and writes model/iJL1143.xml. Every change is a rule applied to the whole
network and is written to ``curation_log.csv`` with its evidence and a
Thiele & Palsson (2010) confidence score (4 biochemical, 3 genetic,
2 sequence/physiological, 1 modelling only, 0 none).

Stages
  1. Strict draft: 1:1 orthologs (tier 1), same-orthogroup genes (tier 2),
     enzyme complexes kept with OLD_ placeholders (RAVEN convention).
  2. Rescue removed template reactions whose EC number is annotated in the
     F. venenatum UniProt proteome (confidence 2), then needs-based gap-fill
     for biomass producibility on minimal glucose medium (confidence 1).
  3. Secretion rule: a metabolite may leave the cell only via a gene-backed
     transporter, the literature/diffusion whitelist, or CO2/H2O/H+.
  4. Redox/energy rules: NADH/NADPH cannot cross the inner mitochondrial
     membrane; mitochondrial ADH runs reductively (Bakker et al. 2000);
     reversible NAD(P) oxidoreductases get eQuilibrator directions
     (|ln Γ| >= 7, Noor et al. 2012). Remaining energy leaks are reported.
  5. F. venenatum biomass (growing-mycelium composition, sources in
     COMPOSITION_SOURCE) and NGAM 1.9 (A. niger
     iMA871) as declared assumptions; biomass molecular weight is reported.
  6. (--rules v2, v3) Thermodynamic directions for all reactions and charged
     organic-acid export; v3 also restricts the mitochondrial GTP/GDP carrier
     to GTP import (Vozza et al. 2004).

Usage (from gem/; the defaults are the inputs and settings of iJL1143):
    python curate_model.py
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import re
import sys
from contextlib import redirect_stdout

import cobra
from cobra.flux_analysis import pfba

import build_draft_model as bdm
from calibrate_model import add_ngam
from egc_test import find_energy_leaks

# Minimal glucose medium.
MINIMAL_GLUCOSE = {
    "BOUNDARY_GLCe": -10.0,
    "BOUNDARY_NH4e": -20.0,
    "BOUNDARY_PIe": -20.0,
    "BOUNDARY_SLFe": -20.0,
    "BOUNDARY_O2e": -20.0,
    "BOUNDARY_H2Oe": -1000.0,
}
ALWAYS_EXPORT = {"BOUNDARY_CO2e", "BOUNDARY_H2Oe", "BOUNDARY_He"}
GROWTH_THRESHOLD = 1e-3
NAD_COFACTORS = {"NAD", "NADH", "NADP", "NADPH", "NADm", "NADHm", "NADPm", "NADPHm",
                 "NADp", "NADHp", "NADPp", "NADPHp"}

# Free cofactor routes checked with all uptakes closed. Each entry is the
# route itself (A -> B); the test maximises its reverse as a dissipation, so a
# positive flux means the network can run the route at no cost.
# expected_zero=True marks uphill routes: with in vivo NAD+/NADH and
# NADPH/NADP+ both high, a free NADH -> NADPH conversion would create
# reducing power for anabolism from nothing, and fungi lack a
# transhydrogenase (Bakker et al. 2001). The downhill NADPH -> NADH direction
# (e.g. the fungal mannitol cycle) and compartment transfers via the
# ethanol-acetaldehyde or malate/oxaloacetate shuttles are reported only.
COFACTOR_ROUTES = [
    ("NADH->NADPH (c)", {"NADH": -1, "NADP": -1, "NAD": 1, "NADPH": 1}, True),
    ("NADPH->NADH (c)", {"NADPH": -1, "NAD": -1, "NADP": 1, "NADH": 1}, False),
    ("NADHm->NADPHm", {"NADHm": -1, "NADPm": -1, "NADm": 1, "NADPHm": 1}, True),
    ("NADPHm->NADHm", {"NADPHm": -1, "NADm": -1, "NADPm": 1, "NADHm": 1}, False),
    ("NADH c->m", {"NADH": -1, "NADm": -1, "NAD": 1, "NADHm": 1}, False),
    ("NADH m->c", {"NADHm": -1, "NAD": -1, "NADm": 1, "NADH": 1}, False),
]
PYRIDINE = {"NAD": "NADP", "NADH": "NADPH"}


# ── Logging ──────────────────────────────────────────────────────────────────

class CurationLog:
    FIELDS = ["reaction_id", "stage", "action", "rule", "evidence", "confidence"]

    def __init__(self) -> None:
        self.rows: list[dict] = []

    def add(self, reaction_id, stage, action, rule, evidence="", confidence=""):
        self.rows.append(dict(reaction_id=reaction_id, stage=stage, action=action,
                              rule=rule, evidence=evidence, confidence=confidence))

    def write(self, path: str) -> None:
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=self.FIELDS)
            w.writeheader()
            w.writerows(self.rows)


# ── Helpers ──────────────────────────────────────────────────────────────────

def set_medium(model: cobra.Model, medium: dict[str, float]) -> None:
    """Close every BOUNDARY_ uptake, then open exactly ``medium``."""
    for r in model.reactions:
        if r.id.startswith("BOUNDARY_"):
            r.lower_bound = 0.0
    for rid, lb in medium.items():
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound = lb


def growth(model: cobra.Model) -> float:
    return model.slim_optimize(error_value=0.0)


def gene_backed(rxn: cobra.Reaction) -> bool:
    """True if the GPR names at least one target-organism gene (not OLD_)."""
    return any(not g.id.startswith("OLD_") for g in rxn.genes)


def is_transport(rxn: cobra.Reaction) -> bool:
    return len({m.compartment for m in rxn.metabolites}) > 1


# ── Stage 1: strict draft ────────────────────────────────────────────────────

def build_strict_draft(template, mapping, orthogroups, log: CurationLog):
    evidence: list[dict] = []
    with redirect_stdout(io.StringIO()):
        draft = bdm.build_draft_model(template, mapping, complex_policy="raven",
                                      fallback=orthogroups, evidence_out=evidence)
    rules = {
        "removed": ("removed", "strict homology: no mapped gene", 0),
        "tier2_orthogroup": ("kept", "tier 2: same OrthoFinder orthogroup", 2),
        "complex_partial": ("kept", "RAVEN: complex kept, unmapped subunits OLD_", 2),
    }
    for e in evidence:
        if e["status"] in rules:
            action, rule, conf = rules[e["status"]]
            log.add(e["reaction_id"], "1_draft", action, rule, e["new_gpr"], conf)
    return draft, evidence


# ── Stage 2: EC rescue and needs-based gap-fill ──────────────────────────────

def load_ec_table(path: str) -> dict[str, list[str]]:
    """UniProt TSV (Entry, EC number, ...) -> {EC: [accessions]} (complete ECs only)."""
    out: dict[str, list[str]] = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            for ec in (row.get("EC number") or "").split(";"):
                ec = ec.strip()
                if re.fullmatch(r"\d+\.\d+\.\d+\.\d+", ec):
                    out.setdefault(ec, []).append(row["Entry"])
    return out


def reaction_ecs(rxn: cobra.Reaction) -> list[str]:
    ec = rxn.annotation.get("ec-code", [])
    ecs = [ec] if isinstance(ec, str) else list(ec)
    return [e for e in ecs if re.fullmatch(r"\d+\.\d+\.\d+\.\d+", e)]


def ec_rescue(model, template, ec_table, log: CurationLog) -> list[str]:
    """Re-add removed template reactions whose every EC is in the proteome."""
    added = []
    for r in template.reactions:
        if r.id in model.reactions:
            continue
        ecs = reaction_ecs(r)
        if not ecs or not all(e in ec_table for e in ecs):
            continue
        genes = sorted({g for e in ecs for g in ec_table[e]})
        new = r.copy()
        model.add_reactions([new])
        model.reactions.get_by_id(r.id).gene_reaction_rule = " or ".join(genes)
        log.add(r.id, "2_ec_rescue", "added", "EC annotated in F. venenatum UniProt proteome",
                f"EC {','.join(ecs)}: {' or '.join(genes)}", 2)
        added.append(r.id)
    return added


def max_atp_yield(model, carbon: str = "BOUNDARY_GLCe") -> float:
    """Max mol ATP per mol carbon source, aerobic minimal medium, no maintenance."""
    with model:
        medium = dict(MINIMAL_GLUCOSE)
        medium.pop("BOUNDARY_GLCe")
        medium.update({carbon: -1.0, "BOUNDARY_O2e": -1000.0})
        set_medium(model, medium)
        for r in model.reactions:
            if r.lower_bound > 0:
                r.lower_bound = 0.0
        atp = cobra.Reaction("ATP_yield_test")
        model.add_reactions([atp])
        atp.add_metabolites({model.metabolites.get_by_id(m): c for m, c in
                             (("ATP", -1), ("H2O", -1), ("ADP", 1), ("PI", 1), ("H", 1))})
        model.objective = atp
        return model.slim_optimize(error_value=0.0)


# Functional tasks a curated network must pass (Thiele & Palsson 2010,
# steps 60-66 biomass precursors and 81-82 known capabilities). The
# respiration threshold sits between substrate-level phosphorylation alone
# (~4-6 ATP/glucose) and the textbook aerobic yield (~30-35): it only asks
# that oxidative phosphorylation works, F. venenatum being a respiring fungus.
RESPIRATION_MIN_ATP = 20.0


def default_tasks(medium):
    def biomass(model):
        with model:
            set_medium(model, medium)
            return growth(model) >= GROWTH_THRESHOLD

    def respiration(model):
        return max_atp_yield(model) >= RESPIRATION_MIN_ATP

    return {"biomass producibility on minimal glucose": biomass,
            f"respiratory ATP (>= {RESPIRATION_MIN_ATP:g} ATP/glucose)": respiration}


def needs_based_gapfill(model, template, tasks, log: CurationLog,
                        exclude: set[str] | frozenset = frozenset(),
                        stage: str = "2_gapfill") -> list[str]:
    """Greedy minimal re-addition of template reactions so all ``tasks`` pass.

    ``tasks`` maps a description to a function(model) -> bool. All missing
    template reactions (minus ``exclude``) are added, then each is removed
    again in template order unless that makes a task fail.
    """
    def all_pass():
        return all(fn(model) for fn in tasks.values())

    if all_pass():
        return []
    cands = [r for r in template.reactions
             if r.id not in model.reactions and r.id not in exclude]
    model.add_reactions([r.copy() for r in cands])
    for r in cands:
        model.reactions.get_by_id(r.id).gene_reaction_rule = ""
    if not all_pass():
        failing = [k for k, fn in tasks.items() if not fn(model)]
        model.remove_reactions([model.reactions.get_by_id(r.id) for r in cands],
                               remove_orphans=True)
        raise RuntimeError(f"Tasks fail even with all template reactions: {failing}")
    kept = []
    for r in cands:
        rx = model.reactions.get_by_id(r.id)
        bounds = rx.bounds
        rx.bounds = (0.0, 0.0)
        if not all_pass():
            rx.bounds = bounds
            kept.append(r.id)
    model.remove_reactions([model.reactions.get_by_id(r.id) for r in cands
                            if r.id not in kept], remove_orphans=True)
    for rid in kept:
        log.add(rid, stage, "added", "needed for functional tasks: " + "; ".join(tasks),
                "template iJB1325, no F. venenatum gene", 1)
    return kept


# ── Stage 3: secretion rule ──────────────────────────────────────────────────

def load_whitelist(path: str) -> dict[str, str]:
    with open(path, newline="") as f:
        return {row["exchange_id"]: f"{row['basis']}; {row['citation']}"
                for row in csv.DictReader(f)}


def apply_secretion_rule(model, whitelist, medium, log: CurationLog,
                         threshold: float = GROWTH_THRESHOLD) -> list[str]:
    """Close exports lacking a gene-backed transporter or literature basis.

    If biomass becomes infeasible, the minimal set of closed exports needed
    for producibility is reopened (confidence 1).
    """
    closed = []
    for ex in [r for r in model.reactions if r.id.startswith("BOUNDARY_")]:
        if ex.upper_bound <= 0 or ex.id in ALWAYS_EXPORT:
            continue
        if ex.id in whitelist:
            log.add(ex.id, "3_secretion", "kept_open", "whitelist", whitelist[ex.id], 2)
            continue
        met = next(iter(ex.metabolites))
        transporters = [r for r in met.reactions if r is not ex and is_transport(r)]
        backed = [r.id for r in transporters if gene_backed(r)]
        if backed:
            log.add(ex.id, "3_secretion", "kept_open", "gene-backed transporter",
                    ",".join(backed), 2)
            continue
        ex.upper_bound = 0.0
        closed.append(ex.id)
        log.add(ex.id, "3_secretion", "export_closed",
                "no gene-backed transporter or literature basis", "", "")
    reopened = []
    with model:
        set_medium(model, medium)
        ok = growth(model) >= threshold
    if not ok:
        for rid in closed:
            model.reactions.get_by_id(rid).upper_bound = 1000.0
        with model:
            set_medium(model, medium)
            for rid in closed:
                rx = model.reactions.get_by_id(rid)
                rx.upper_bound = 0.0
                if growth(model) < threshold:
                    rx.upper_bound = 1000.0
                    reopened.append(rid)
        for rid in closed:
            model.reactions.get_by_id(rid).upper_bound = 1000.0 if rid in reopened else 0.0
        for rid in reopened:
            log.add(rid, "3_secretion", "export_reopened",
                    "needed for biomass producibility", "", 1)
    return [c for c in closed if c not in reopened]


# ── Stage 4: redox / energy rules ────────────────────────────────────────────

def block_pyridine_nucleotide_transport(model, log: CurationLog) -> list[str]:
    """NADH/NADPH do not cross the inner mitochondrial membrane (Bakker 2001)."""
    blocked = []
    for r in model.reactions:
        ids = {m.id for m in r.metabolites}
        if (is_transport(r) and not gene_backed(r)
                and any(i in ids for i in ("NADH", "NADPH"))
                and any(i in ids for i in ("NADHm", "NADPHm"))):
            r.bounds = (0.0, 0.0)
            blocked.append(r.id)
            log.add(r.id, "4_redox", "blocked",
                    "NAD(P)H cannot cross the inner mitochondrial membrane",
                    "Bakker et al. 2001 FEMS Microbiol Rev 25:15-37", 3)
    return blocked


def mitochondrial_adh_reductive(model, log: CurationLog) -> list[str]:
    """Mitochondrial ADH reoxidises matrix NADH (acetaldehyde -> ethanol),
    i.e. the ethanol-acetaldehyde shuttle runs mitochondria -> cytosol
    (Bakker et al. 2000 J Bacteriol 182:4730)."""
    changed = []
    for r in model.reactions:
        if "1.1.1.1" not in reaction_ecs(r):
            continue
        ids = {m.id: c for m, c in r.metabolites.items()}
        if not {"ETHm", "ACALm", "NADm", "NADHm"} <= set(ids):
            continue
        # forward produces NADHm if NADHm has a positive coefficient
        if ids["NADHm"] > 0:
            r.upper_bound = min(r.upper_bound, 0.0)
        else:
            r.lower_bound = max(r.lower_bound, 0.0)
        changed.append(r.id)
        log.add(r.id, "4_redox", "direction_set", "mitochondrial ADH reductive only",
                "Bakker et al. 2000 J Bacteriol 182:4730-4737", 3)
    return changed


def _eq_compound_id(met: cobra.Metabolite) -> str | None:
    kegg = met.annotation.get("kegg.compound")
    if kegg:
        return "kegg:" + (kegg if isinstance(kegg, str) else kegg[0])
    chebi = met.annotation.get("chebi")
    if chebi:
        return "chebi:" + (chebi if isinstance(chebi, str) else chebi[0])
    return None


def thermodynamic_directions(model, medium, log: CurationLog, cc=None,
                             ln_gamma_min: float = 7.0, redox_only: bool = True,
                             stage: str = "4_thermo") -> dict[str, int]:
    """Assign directions to reversible single-compartment NAD(P) oxidoreductases
    (``redox_only``, the v1 rule) or to every reversible metabolic reaction
    (v2: the same test, no longer limited to NAD(P) reactions).

    A reaction is made irreversible when |ln Γ| >= ``ln_gamma_min`` and the
    estimate is more than two standard errors from zero (Noor et al. 2012).
    A direction that would make biomass infeasible is reverted and flagged.
    Returns counts per outcome.
    """
    if cc is None:
        from equilibrator_api import ComponentContribution
        cc = ComponentContribution()
    counts = {"evaluated": 0, "set": 0, "reversible": 0, "unmapped": 0, "conflict": 0}
    for r in model.reactions:
        if r.lower_bound >= 0 or r.upper_bound <= 0 or r.boundary or is_transport(r):
            continue
        if redox_only and not ({m.id for m in r.metabolites} & NAD_COFACTORS):
            continue
        parts = []
        for m, coeff in r.metabolites.items():
            if m.id in ("H", "Hm", "Hp"):
                continue
            cid = _eq_compound_id(m)
            if cid is None:
                parts = None
                break
            parts.append((coeff, cid))
        if not parts:
            counts["unmapped"] += 1
            continue
        lhs = " + ".join(f"{-c:g} {cid}" for c, cid in parts if c < 0)
        rhs = " + ".join(f"{c:g} {cid}" for c, cid in parts if c > 0)
        try:
            rxn = cc.parse_reaction_formula(f"{lhs} = {rhs}")
            if not rxn.is_balanced(ignore_atoms=("H",)):
                counts["unmapped"] += 1
                continue
            lng = cc.ln_reversibility_index(rxn)
            val, err = float(lng.value.magnitude), float(lng.error.magnitude)
        except Exception:
            counts["unmapped"] += 1
            continue
        counts["evaluated"] += 1
        if abs(val) < ln_gamma_min or abs(val) < 2 * err:
            counts["reversible"] += 1
            continue
        old = r.bounds
        if val > 0:   # forward strongly unfavourable -> reverse only
            r.upper_bound = 0.0
            direction = "reverse"
        else:
            r.lower_bound = 0.0
            direction = "forward"
        with model:
            set_medium(model, medium)
            feasible = growth(model) >= GROWTH_THRESHOLD
        if not feasible:
            r.bounds = old
            counts["conflict"] += 1
            log.add(r.id, stage, "conflict_kept_reversible",
                    "eQuilibrator direction would block biomass",
                    f"ln_Gamma={val:.1f}+/-{err:.1f}", "")
            continue
        counts["set"] += 1
        log.add(r.id, stage, f"direction_{direction}",
                "eQuilibrator |ln Gamma|>=7 (Noor 2012)",
                f"ln_Gamma={val:.1f}+/-{err:.1f}", 2)
    return counts


# Net negative charge at pH 7 of organic acids that have a cytosol <-> extracellular
# transport in iJB1325, from their pKa values (Dawson et al., Data for Biochemical
# Research, 3rd ed.). Amino acids count by their net charge (Asp, Glu: -1).
ACID_CHARGE_PH7 = {
    "AC": 1, "PYR": 1, "LAC": 1, "LLAC": 1, "FOR": 1, "PROP": 1, "GLCNT": 1,
    "GALUNT": 1, "GLCUNT": 1, "URIC": 1, "ATT": 1, "ASP": 1, "GLU": 1,
    "SUCC": 2, "FUM": 2, "MAL": 2, "AKG": 2, "OA": 2, "OXAL": 2, "TAR": 2,
    "CIT": 3, "ICIT": 3, "ACO": 3,
}
ATP_PER_PROTON = 1.0   # plasma-membrane H+-ATPase, 1 H+/ATP (Serrano 1988 BBA 947:1)


def restrict_gtp_carrier(model, log: CurationLog) -> None:
    """The mitochondrial GTP/GDP carrier imports GTP in exchange for GDP (Vozza et al. 2004, Ggc1p).
    Open both ways it lets mitochondrial ATP reach the cytosol as GTP (adenylate kinases r541, r542),
    without the proton that the ADP/ATP carrier r251 costs, which sets P/O above the template's design."""
    rxn = model.reactions.get_by_id("r2335")          # GDPm + GTP <=> GDP + GTPm
    rxn.bounds = (0.0, rxn.upper_bound)
    log.add("r2335", "6_v3_gtp_carrier", "direction_set", "GTP/GDP carrier imports GTP",
            "Vozza et al. 2004, J Biol Chem 279:20850", 3)


def charge_acid_export(model, charges: dict[str, int], log: CurationLog,
                       atp_per_h: float = ATP_PER_PROTON) -> list[str]:
    """Make organic-acid export pay for the protons it leaves behind.

    The model carries acids as neutral molecules, so exporting one is free.
    In the cell the anion leaves and its protons stay in the cytosol, where
    the plasma-membrane H+-ATPase must pump them out (van Maris et al. 2004
    Metab Eng 6:245). Each cytosol <-> extracellular transport of an acid in
    ``charges`` (id -> net negative charge at pH 7) keeps its inward
    direction, which stays free (the protonated acid diffuses in), and gains
    a separate export reaction ``<id>_exp`` consuming charge x atp_per_h ATP.
    Returns the ids of the transports changed.
    """
    atp, adp, pi, h2o = (model.metabolites.get_by_id(m) for m in ("ATP", "ADP", "PI", "H2O"))
    changed = []
    for cyt_id, n in charges.items():
        if cyt_id not in model.metabolites or f"{cyt_id}e" not in model.metabolites:
            continue
        cyt, ext = model.metabolites.get_by_id(cyt_id), model.metabolites.get_by_id(f"{cyt_id}e")
        for t in list(cyt.reactions & ext.reactions):
            if set(t.metabolites) != {cyt, ext}:
                continue
            if t.metabolites[cyt] > 0:          # written inward: keep the forward direction
                t.lower_bound = max(t.lower_bound, 0.0)
            else:                               # written outward: keep the reverse direction
                t.upper_bound = min(t.upper_bound, 0.0)
            exp = cobra.Reaction(f"{t.id}_exp", name=f"{cyt.name} export with H+ extrusion",
                                 lower_bound=0.0, upper_bound=1000.0)
            model.add_reactions([exp])
            cost = n * atp_per_h
            exp.add_metabolites({cyt: -1, ext: 1, atp: -cost, h2o: -cost, adp: cost, pi: cost})
            log.add(t.id, "6_v2_acid_export", "export_charged",
                    "anion export: protons extruded by PM H+-ATPase (van Maris 2004)",
                    f"{n} H+ per {cyt_id}, {atp_per_h:g} ATP per H+ -> {exp.id}", 2)
            changed.append(t.id)
    return changed


def remove_cofactor_duplicates(model, log: CurationLog) -> list[str]:
    """Remove a gene-less NAD(H) or NADP(H) variant of a gene-backed reaction.

    Two reactions that differ only by NAD(H) <-> NADP(H) form a free
    transhydrogenase cycle. When only one of them has F. venenatum gene
    evidence, the gene-less variant is removed (cf. Agren et al. 2013,
    iAL1006, which blocked NADH/NADPH futile cycles).
    """
    def key(r, swap):
        out = []
        for m, c in r.metabolites.items():
            mid = m.id
            for nad, nadp in PYRIDINE.items():
                for suffix in ("", "m", "p"):
                    if swap and mid == nadp + suffix:
                        mid = nad + suffix
            out.append((mid, round(c, 6)))
        return frozenset(out)

    by_key: dict[frozenset, list] = {}
    for r in model.reactions:
        if r.boundary:
            continue
        by_key.setdefault(key(r, swap=True), []).append(r)
    removed = []
    for group in by_key.values():
        if len(group) < 2:
            continue
        backed = [r for r in group if gene_backed(r)]
        geneless = [r for r in group if not r.genes]
        if not backed or not geneless:
            continue
        for r in geneless:
            if key(r, swap=False) in {key(b, swap=False) for b in backed}:
                continue  # identical isoenzyme, not a cofactor variant
            removed.append(r.id)
            log.add(r.id, "4_redox", "removed",
                    "gene-less cofactor variant of a gene-backed reaction",
                    "gene-backed variant: " + ",".join(b.id for b in backed), 2)
    model.remove_reactions([model.reactions.get_by_id(i) for i in removed])
    return removed


def orient_cofactor_pairs(model, medium, log: CurationLog) -> list[str]:
    """Orient NAD/NADP isoenzyme pairs by in vivo cofactor ratios.

    Reactions that differ only by NAD(H) <-> NADP(H) and both have gene
    evidence still form a free transhydrogenase cycle. In vivo NAD+/NADH is
    high and NADPH/NADP+ is high, so NAD-linked variants run oxidatively and
    NADP-linked variants reductively (Pereira, Nielsen & Rocha 2016, Metab Eng
    Commun 3:153). Applied to every such pair; an orientation that would make
    biomass infeasible is reverted and flagged.
    """
    groups: dict[frozenset, list] = {}
    for r in model.reactions:
        if r.boundary:
            continue
        key = []
        for m, c in r.metabolites.items():
            mid = m.id
            for nad, nadp in PYRIDINE.items():
                for suffix in ("", "m", "p"):
                    if mid == nadp + suffix:
                        mid = nad + suffix
            key.append((mid, round(c, 6)))
        groups.setdefault(frozenset(key), []).append(r)
    with model:
        set_medium(model, medium)
        base_ok = growth(model) >= GROWTH_THRESHOLD
    changed = []
    for group in groups.values():
        ids = [{m.id for m in r.metabolites} for r in group]
        has_nad = any(i & {"NADH", "NADHm", "NADHp"} for i in ids)
        has_nadp = any(i & {"NADPH", "NADPHm", "NADPHp"} for i in ids)
        if len(group) < 2 or not (has_nad and has_nadp):
            continue
        # only fully reversible pairs form a free cycle; irreversible
        # reactions keep their curated direction
        if not all(r.lower_bound < 0 < r.upper_bound for r in group):
            continue
        for r in group:
            red = next(m for m in r.metabolites
                       if m.id in ("NADH", "NADHm", "NADHp", "NADPH", "NADPHm", "NADPHp"))
            nadp_linked = "NADP" in red.id
            produces_red_forward = r.metabolites[red] > 0
            old = r.bounds
            # NADP-linked: consume NADPH (reductive); NAD-linked: produce NADH (oxidative)
            if nadp_linked == produces_red_forward:
                r.lower_bound, r.upper_bound = min(r.lower_bound, 0.0), 0.0
                direction = "reverse"
            else:
                r.lower_bound, r.upper_bound = 0.0, max(r.upper_bound, 0.0)
                direction = "forward"
            with model:
                set_medium(model, medium)
                ok = growth(model) >= GROWTH_THRESHOLD
            if base_ok and not ok:
                r.bounds = old
                log.add(r.id, "4_redox", "conflict_kept_reversible",
                        "cofactor-ratio orientation would block biomass", "", "")
                continue
            changed.append(r.id)
            log.add(r.id, "4_redox", f"direction_{direction}",
                    ("NADP-linked runs reductively" if nadp_linked else
                     "NAD-linked runs oxidatively") + " (in vivo cofactor ratios)",
                    "Pereira et al. 2016 Metab Eng Commun 3:153-163; pair: "
                    + ",".join(x.id for x in group), 2)
    return changed


def remove_unsupported_transhydrogenases(model, ec_table, log: CurationLog) -> list[str]:
    """Remove NAD(P) transhydrogenation unless a gene is annotated EC 1.6.1.-.

    Fungi lack pyridine nucleotide transhydrogenase (Anderlund et al. 1999;
    Bakker et al. 2001), so a template GPR alone is not enough evidence.
    """
    ec_genes = {g for ec, genes in ec_table.items() if ec.startswith("1.6.1.") for g in genes}
    removed = []
    def base_id(mid: str) -> str:
        if mid[-1:] in ("m", "p") and mid[:-1] in ("NAD", "NADH", "NADP", "NADPH", "H"):
            return mid[:-1]
        return mid

    for r in list(model.reactions):
        base = {base_id(m.id) for m in r.metabolites} - {"H"}
        if base != {"NAD", "NADH", "NADP", "NADPH"}:
            continue
        if any(g.id in ec_genes for g in r.genes):
            continue
        removed.append(r.id)
        log.add(r.id, "4_redox", "removed",
                "transhydrogenase without EC 1.6.1.- annotated gene",
                f"GPR was '{r.gene_reaction_rule}'; Anderlund 1999, Bakker 2001", 2)
    model.remove_reactions([model.reactions.get_by_id(i) for i in removed])
    return removed


def report_leaks(model) -> list[tuple[str, float, bool]]:
    """Fritzemeier EGC tests plus free-cofactor-route tests.

    Returns (label, max_flux, expected_zero).
    """
    tests = [(f"route {label}", {k: -v for k, v in st.items()})
             for label, st, _ in COFACTOR_ROUTES]
    expected = {f"route {label}": exp for label, _, exp in COFACTOR_ROUTES}
    return [(label, flux, expected.get(label, True))
            for label, flux in find_energy_leaks(model, extra_tests=tests)]


# ── Stage 5: biomass and maintenance ─────────────────────────────────────────

# Generic lipid formulas carry R groups. R is approximated by a C18:1 acyl
# chain (C17H33): this reproduces dioleoyl-PE (744 g/mol) and DOPC (786 g/mol).
R_GROUP = {"C": 17, "H": 33}
_ATOMIC_WEIGHT = {"C": 12.011, "H": 1.008, "O": 15.999, "N": 14.007, "P": 30.974,
                  "S": 32.06, "Fe": 55.845, "Mg": 24.305, "K": 39.098, "Ca": 40.078,
                  "Na": 22.990, "Cl": 35.45, "Co": 58.933, "Cu": 63.546, "Zn": 65.38,
                  "Mn": 54.938, "Se": 78.971, "Mo": 95.95}


def formula_weight(met: cobra.Metabolite, r_group=R_GROUP) -> float | None:
    """g/mol from the formula, with R groups replaced by ``r_group``."""
    if not met.formula:
        return None
    elements = dict(met.elements)
    n_r = elements.pop("R", 0)
    for el, n in r_group.items():
        elements[el] = elements.get(el, 0) + n * n_r
    if any(el not in _ATOMIC_WEIGHT for el in elements):
        return None
    return sum(_ATOMIC_WEIGHT[el] * n for el, n in elements.items())


def biomass_weight(model, biomass_met: str = "BIOMASS") -> float | None:
    """g/mmol of BIOMASS, expanding formula-less pseudo-metabolites through
    their single producing reaction. None if not computable."""
    cache: dict[str, float | None] = {}

    def weight(met: cobra.Metabolite, depth: int = 0) -> float | None:
        if met.id in cache:
            return cache[met.id]
        fw = formula_weight(met)
        if fw is not None:
            w = fw / 1000.0  # g/mmol
        elif depth > 5:
            w = None
        else:
            prod = [r for r in met.reactions if r.metabolites[met] > 0 and not r.boundary]
            if len(prod) != 1:
                w = None
            else:
                r, total = prod[0], 0.0
                for m, c in r.metabolites.items():
                    if m is met:
                        continue
                    mw = weight(m, depth + 1)
                    if mw is None:
                        total = None
                        break
                    total -= c * mw
                w = None if total is None else total / r.metabolites[met]
        cache[met.id] = w
        return w

    if biomass_met not in model.metabolites:
        return None
    return weight(model.metabolites.get_by_id(biomass_met))


# F. venenatum A3/5 biomass composition (% of ash-free dry mass) for growing mycelium (version 1.1.0).
# Measured on dry mycelium: protein 37 % and RNA 9 % (US patent 4,041,189: 37.3 % protein and 8.58 % RNA in the
# starting mycelium, RNA 8.7 +/- 1.1 % in untreated controls; Trinci 1992: RNA 8-9 % at 0.19 h-1).
# Cell wall 17 % and lipid 9 %: the mycoprotein product values (FDA response letter to GRN 000091, 2002:
# cell-wall material 22-28 %, fat 12-14 % of dry weight) times the 0.68 g of product left per g of mycelium after
# RNA reduction (US 4,041,189), assuming walls and lipids are retained. DNA 0.9 %: A. oryzae chemostat
# (Pedersen et al. 1999). Ash 5 % (US 3,937,654; product specification <= 5 %) is removed and the values are
# expressed on ash-free mass. The small-molecule pool is the remainder.
BIOMASS_COMPOSITION = {"Protein": 38.9, "CELLWALL": 17.9, "LIPIDS": 9.5, "RNA": 9.5,
                       "DNA": 0.9, "POOL": 23.3}
COMPOSITION_SOURCE = {
    "Protein": "37.3% of dry mycelium (US 4041189); ash-free",
    "CELLWALL": "product cell-wall material 22-28% (FDA GRN 000091 letter 2002) x 0.68 retained (US 4041189); ash-free",
    "LIPIDS": "product fat 12-14% (FDA GRN 000091 letter 2002) x 0.68 retained (US 4041189); ash-free",
    "RNA": "8.7% in untreated mycelium (US 4041189); 8-9% (Trinci 1992); ash-free",
    "DNA": "A. oryzae chemostat (Pedersen et al. 1999)",
    "POOL": "remainder of ash-free dry mass",
}


def normalize_biomass(model, template, fractions, log: CurationLog,
                      growth_id: str = "GROWTH") -> dict[str, float]:
    """Set GROWTH coefficients so each macromolecule has its target mass
    fraction and BIOMASS weighs 1 g/mmol (Chan et al. 2017).

    Per-mmol weights come from formulas, expanding pseudo-metabolites through
    their assembly reactions. A component whose weight cannot be computed
    is inferred from the template, assuming the template GROWTH reaction is
    itself normalised to 1 g/mmol (not needed once R groups are approximated).
    """
    total = sum(fractions.values())
    g_model = model.reactions.get_by_id(growth_id)
    g_tpl = template.reactions.get_by_id(growth_id)
    weights = {c: biomass_weight(model, c) for c in fractions}
    missing = [c for c, w in weights.items() if w is None]
    if len(missing) > 1:
        raise ValueError(f"Cannot infer weights of several components: {missing}")
    for c in missing:
        known = sum(-coeff * biomass_weight(template, m.id)
                    for m, coeff in g_tpl.metabolites.items()
                    if coeff < 0 and m.id in fractions and m.id != c)
        tpl_coeff = -g_tpl.metabolites[template.metabolites.get_by_id(c)]
        weights[c] = (1.0 - known) / tpl_coeff
        log.add(growth_id, "5_biomass", "inferred_weight",
                f"{c} g/mmol inferred assuming template biomass = 1 g/mmol",
                f"{weights[c]:.4f}", 1)
    new = {}
    for c, frac in fractions.items():
        met = model.metabolites.get_by_id(c)
        coeff = (frac / total) / weights[c]
        g_model.add_metabolites({met: -coeff - g_model.metabolites[met]})
        new[c] = coeff
        log.add(growth_id, "5_biomass", "coefficient_set",
                f"{c}: {frac / total:.3f} g/gDW at {weights[c]:.4f} g/mmol",
                f"{COMPOSITION_SOURCE.get(c, 'see curate_model.py')} (renormalised from {total:.1f}%) "
                f"-> {coeff:.4g} mmol", 2)
    return new


# ── Main ─────────────────────────────────────────────────────────────────────

def summarize(model, label: str, medium) -> str:
    with model:
        set_medium(model, medium)
        mu = growth(model)
        ex_open = sum(1 for r in model.reactions
                      if r.id.startswith("BOUNDARY_") and r.upper_bound > 0)
    n_gpr = sum(1 for r in model.reactions if r.gene_reaction_rule)
    return (f"{label}: {len(model.reactions)} rxns, {len(model.genes)} genes, "
            f"{n_gpr} with GPR, {ex_open} open exports, mu(glc 10)={mu:.4f}")


HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    p.add_argument("--template", default=os.path.join(HERE, "raw_data", "MOESM3_cbs51388.xml"))
    p.add_argument("--orthologs", default=os.path.join(HERE, "inputs", "orthologs_cbs513.tsv"))
    p.add_argument("--orthogroups", default=os.path.join(HERE, "inputs", "Orthogroups.tsv"))
    p.add_argument("--ec-table", default=os.path.join(HERE, "inputs", "fven_uniprot_ec.tsv"),
                   help="UniProt TSV for proteome UP000245910 (Entry, EC number, ...)")
    p.add_argument("--whitelist", default=os.path.join(HERE, "inputs", "secretion_whitelist.csv"))
    p.add_argument("--ngam", type=float, default=1.9)
    p.add_argument("--template-species", default="Aniger_CBS513",
                   help="column of the template genes in Orthogroups.tsv")
    p.add_argument("--rules", choices=("v1", "v2", "v3"), default="v3",
                   help="v2 adds stage 6: thermodynamic directions for every "
                        "reaction and charged organic-acid export; v3 also "
                        "restricts the mitochondrial GTP/GDP carrier to GTP import")
    p.add_argument("--output", default=os.path.join(HERE, "build", "curated.xml"))
    p.add_argument("--log-dir", default=os.path.join(HERE, "..", "model", "curation"))
    args = p.parse_args()

    os.makedirs(args.log_dir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    log = CurationLog()
    template = cobra.io.read_sbml_model(args.template)
    medium = MINIMAL_GLUCOSE

    # 1
    model, evidence = build_strict_draft(
        template, bdm.load_ortholog_mapping(args.orthologs),
        bdm.load_orthogroups(args.orthogroups, template_col=args.template_species), log)
    with open(os.path.join(args.log_dir, "draft_evidence.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(evidence[0]), delimiter="\t")
        w.writeheader()
        w.writerows(evidence)
    print(summarize(model, "1 strict draft", medium))

    # 2
    ec_table = load_ec_table(args.ec_table)
    rescued = ec_rescue(model, template, ec_table, log)
    print(f"  EC rescue added {len(rescued)}")
    tasks = default_tasks(medium)
    filled = needs_based_gapfill(model, template, tasks, log)
    print(f"  gap-fill added {len(filled)}")
    with template:
        set_medium(template, medium)
        mu_tpl = growth(template)
    print(summarize(model, "2 rescued + gap-filled", medium) + f"  (template mu={mu_tpl:.4f})")

    # 3
    closed = apply_secretion_rule(model, load_whitelist(args.whitelist), medium, log)
    print(f"  exports closed {len(closed)}")
    print(summarize(model, "3 secretion rule", medium))

    # 4
    nb = block_pyridine_nucleotide_transport(model, log)
    na = mitochondrial_adh_reductive(model, log)
    nd = remove_cofactor_duplicates(model, log)
    nt = remove_unsupported_transhydrogenases(model, ec_table, log)
    no = orient_cofactor_pairs(model, medium, log)
    counts = thermodynamic_directions(model, medium, log)
    removed_by_rule = {r["reaction_id"] for r in log.rows
                       if r["stage"] == "4_redox" and r["action"] in ("removed", "blocked")}
    refill = needs_based_gapfill(model, template, tasks, log,
                                 exclude=removed_by_rule, stage="4_task_recheck")
    print(f"  task re-check after redox rules added {refill}; "
          f"max ATP/glucose = {max_atp_yield(model):.1f}")
    print(f"  NAD(P)H transport blocked {nb}, mito ADH set {na}, cofactor variants "
          f"removed {nd}, transhydrogenases removed {nt}, "
          f"cofactor pairs oriented {no}, thermo {counts}")
    print(summarize(model, "4 redox rules", medium))
    leaks = report_leaks(model)
    with open(os.path.join(args.log_dir, "energy_leaks.tsv"), "w") as f:
        f.write("test\tmax_flux\texpected_zero\n")
        for label, flux, exp in leaks:
            f.write(f"{label}\t{flux:.6g}\t{exp}\n")
    for label, flux, exp in leaks:
        if flux == flux and abs(flux) > 1e-6:
            print(f"  {'LEAK' if exp else 'allowed route (reported)'}: {label} = {flux:.4g}")

    # 5
    normalize_biomass(model, template, BIOMASS_COMPOSITION, log)
    bw = biomass_weight(model)
    log.add("GROWTH", "5_biomass", "checked",
            "biomass molecular weight after normalisation (g/mmol); lipid R groups "
            "approximated as C18:1 acyl (C17H33)", "n/a" if bw is None else f"{bw:.3f}", "")
    print(f"  biomass normalised to the growing-mycelium composition (1 g/mmol by construction)")
    with redirect_stdout(io.StringIO()):
        add_ngam(model, args.ngam)
    log.add("NGAM", "5_biomass", "set", f"NGAM={args.ngam} mmol ATP/gDW/h",
            "A. niger iMA871 (Andersen 2008) - assumption", 1)
    print(summarize(model, "5 biomass + NGAM", medium))

    if args.rules in ("v2", "v3"):
        # 6
        thermo2 = thermodynamic_directions(model, medium, log, redox_only=False,
                                           stage="6_v2_thermo")
        charged = charge_acid_export(model, ACID_CHARGE_PH7, log)
        print(f"  v2 thermo (all reactions) {thermo2}; acid exports charged {len(charged)}; "
              f"max ATP/glucose = {max_atp_yield(model):.1f}")
        if args.rules == "v3":
            restrict_gtp_carrier(model, log)
            print(f"  v3 GTP/GDP carrier imports GTP; max ATP/glucose = {max_atp_yield(model):.1f}")
        print(summarize(model, "6 rules", medium))

    set_medium(model, medium)
    model.id = "fven_a3_5_curated"
    model.name = "Fusarium venenatum A3/5 GEM, strict homology + rule-based curation"
    cobra.io.write_sbml_model(model, args.output)
    log.write(os.path.join(args.log_dir, "curation_log.csv"))
    print(f"\nSaved {args.output}; {len(log.rows)} log rows in {args.log_dir}/curation_log.csv")


if __name__ == "__main__":
    main()
