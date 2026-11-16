"""Side-by-side checks of the Draft and Curated (iJL1143) models, before any data comparison.

  draft    iJB1325 (CBS 513.88) transferred unchanged, same biomass and NGAM, no curation
  curated  iJL1143: stages 1-6 of the curation rules

Checks: size; GPR evidence (does a GPR gene's UniProt EC match the reaction
EC; named central-carbon reactions); aerobic biomass yield per sugar; maximum
ATP per glucose; the energy-leak tests of curate_model; the O2-limited
byproduct profile on glucose.

Outputs: results/model_versions_qc/{qc.json, provenance.json}.

Usage (from analysis/):  python scripts/model_versions_qc.py
"""

from __future__ import annotations

import csv
import json
import re

import curate_model as cm
from fbaval import simulate
from fbaval.provenance import write_provenance

from fbaval.paths import CURATED, DRAFT, GEM_INPUTS, RESULTS

OUT = RESULTS / "model_versions_qc"
MODELS = {"draft": DRAFT, "curated": CURATED}
CENTRAL = {"r5a": "hexokinase", "r12a": "FBP aldolase", "r14": "GAPDH", "r15": "phosphoglycerate kinase",
           "r17": "enolase", "r18": "pyruvate kinase", "r31": "citrate synthase",
           "r22": "G6P dehydrogenase", "r1273": "phosphoketolase", "r60": "acetate kinase"}
SUGARS = ("glucose", "fructose", "sucrose", "xylose")


def gene_ecs() -> tuple[dict, dict]:
    ecs, names = {}, {}
    with open(GEM_INPUTS / "fven_uniprot_ec.tsv", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            ecs[r["Entry"]] = {e.strip() for e in (r["EC number"] or "").split(";")
                               if re.fullmatch(r"\d+\.\d+\.\d+\.\d+", e.strip())}
            names[r["Entry"]] = r["Protein names"].split("(")[0].strip()[:60]
    return ecs, names


def gpr_evidence(model, ecs, names) -> dict:
    ok = conflict = unannotated = 0
    for r in model.reactions:
        rxn_ec = set(cm.reaction_ecs(r))
        if not rxn_ec or not r.genes:
            continue
        g_ec = set().union(*(ecs.get(g.id, set()) for g in r.genes))
        if not g_ec:
            unannotated += 1
        elif g_ec & rxn_ec:
            ok += 1
        else:
            conflict += 1
    central = {}
    for rid, label in CENTRAL.items():
        if rid in model.reactions:
            r = model.reactions.get_by_id(rid)
            central[label] = {"bounds": list(r.bounds),
                              "genes": [names.get(g.id, g.id) for g in r.genes][:3]}
    return {"ec_consistent": ok, "ec_conflict": conflict, "gene_ec_unannotated": unannotated,
            "central": central}


def yields(model) -> dict:
    out = {}
    bw = simulate.biomass_weight(model)
    for s in SUGARS:
        with model:
            simulate.set_sugars(model, {s: 1.0})
            mu = model.slim_optimize(error_value=0.0) * bw
        out[s] = mu / (simulate.MW[s] / 1000.0)          # g biomass / g sugar
    return out


def o2_profile(model) -> list[dict]:
    df = simulate.static_profile(model, "glucose", 1.0, o2_fractions=(1.0, 0.6, 0.3))
    keep = ["o2_fraction", "yield_gg"] + [c for c in df if c.endswith(("_min", "_max"))]
    return df[keep].round(4).to_dict(orient="records")


def main() -> None:
    ecs, names = gene_ecs()
    qc = {}
    for label, fname in MODELS.items():
        m = simulate.load(str(fname))
        with m:
            cm.set_medium(m, cm.MINIMAL_GLUCOSE)
            atp = cm.max_atp_yield(m)
        leaks = [(t, f) for t, f, exp in cm.report_leaks(m) if exp and f == f and abs(f) > 1e-6]
        qc[label] = {"file": fname.name, "reactions": len(m.reactions), "genes": len(m.genes),
                     "open_exports": sum(1 for r in m.boundary if r.upper_bound > 0),
                     "biomass_weight": simulate.biomass_weight(m),
                     "gpr": gpr_evidence(m, ecs, names), "yield_gg": yields(m),
                     "max_atp_per_glucose": atp, "energy_leaks": leaks, "o2_profile": o2_profile(m)}
        g = qc[label]["gpr"]
        print(f"{label}: {len(m.reactions)} rxns, {len(m.genes)} genes | EC ok {g['ec_consistent']} "
              f"conflict {g['ec_conflict']} | yield glc {qc[label]['yield_gg']['glucose']:.3f} | "
              f"ATP/glc {atp:.1f} | leaks {leaks}")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "qc.json").write_text(json.dumps(qc, indent=1))
    write_provenance(OUT, params={"models": {k: v.name for k, v in MODELS.items()}, "fitted_parameters": "none"})


if __name__ == "__main__":
    main()
