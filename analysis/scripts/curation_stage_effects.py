"""Which curation stage changes energetics and yield. Two ablations from the curation change log.

Curated without stage k: every reaction the log lists under stage k goes back to its Draft state.
Draft with stage k only: every reaction the log lists under stage k takes its Curated state.
Measures: maximum ATP per glucose, failed energy-generating-cycle tests, open exports and aerobic yield on
glucose and fructose at 1 mmol/gDW/h. Stage 5 (biomass, maintenance) is shared by both models, so it is not ablated.

Outputs: results/curation_stage_effects/{effects.csv, provenance.json}

Usage (from analysis/):  python scripts/curation_stage_effects.py
"""

from __future__ import annotations

import curate_model as cm
import pandas as pd

from fbaval import simulate
from fbaval.provenance import write_provenance

from fbaval.paths import CURATED, CURATION_DIR, DRAFT, RESULTS

OUT = RESULTS / "curation_stage_effects"
ACTIONS_SKIPPED = {"kept", "kept_open", "checked", "coefficient_set", "set", "inferred_weight"}
STAGES = {"1_draft": "1 Strict-homology draft", "2_ec_rescue": "2 EC rescue", "2_gapfill": "2 Gap-fill", "3_secretion": "3 Secretion rule",
          "4_redox": "4 Redox rules", "4_thermo": "4 NAD(P) thermodynamic directions", "6_v2_thermo": "6 Thermodynamic directions, all reactions",
          "6_v2_acid_export": "6 Organic-acid export cost", "6_v3_gtp_carrier": "6 GTP/GDP carrier direction"}


def transfer(target, source, rids) -> dict:
    """Give the listed reactions in ``target`` their state in ``source``: stoichiometry and bounds, adding or removing as needed."""
    n = {"changed": 0, "added": 0, "removed": 0}
    for rid in rids:
        in_t, in_s = rid in target.reactions, rid in source.reactions
        if in_s and in_t:
            rt, rs = target.reactions.get_by_id(rid), source.reactions.get_by_id(rid)
            rt.subtract_metabolites(rt.metabolites)
            rt.add_metabolites({(target.metabolites.get_by_id(m.id) if m.id in target.metabolites else m.copy()): c for m, c in rs.metabolites.items()})
            rt.bounds = rs.bounds
            n["changed"] += 1
        elif in_s:
            r = source.reactions.get_by_id(rid).copy()
            target.add_reactions([r])
            n["added"] += 1
        elif in_t:
            target.remove_reactions([target.reactions.get_by_id(rid)], remove_orphans=False)
            n["removed"] += 1
    return n


def measure(m) -> dict:
    with m:
        cm.set_medium(m, cm.MINIMAL_GLUCOSE)
        atp = cm.max_atp_yield(m)
    leaks = [t for t, f, exp in cm.report_leaks(m) if exp and f == f and abs(f) > 1e-6]
    out = {"max_atp_per_glucose": atp, "cycle_tests_failed": len(leaks), "cycle_tests": "; ".join(leaks),
           "open_exports": sum(1 for r in m.boundary if r.upper_bound > 0)}
    for s in ("glucose", "fructose"):
        with m:
            simulate.set_sugars(m, {s: 1.0})
            mu = m.slim_optimize(error_value=0.0) * simulate.biomass_weight(m)
        out[f"yield_{s}_gg"] = mu / (simulate.MW[s] / 1000.0)
    return out


def main() -> None:
    draft = simulate.load(str(DRAFT), normalize=True)
    cur = simulate.load(str(CURATED), normalize=True)
    log = pd.read_csv(CURATION_DIR / "curation_log.csv")
    log = log[~log.action.isin(ACTIONS_SKIPPED)]
    rows = [{"variant": "Draft", "stage": "", **measure(draft)}, {"variant": "Curated", "stage": "", **measure(cur)}]
    for stage, label in STAGES.items():
        rids = sorted(set(log[log.stage == stage].reaction_id))
        m = cur.copy()
        n = transfer(m, draft, rids)
        rows.append({"variant": "Curated without stage", "stage": label, "n_reactions": len(rids), **n, **measure(m)})
        m = draft.copy()
        n = transfer(m, cur, rids)
        rows.append({"variant": "Draft with stage only", "stage": label, "n_reactions": len(rids), **n, **measure(m)})
    # Combinations: what the Curated ceiling depends on, and what brings the Draft down to it
    for variant, label, stages, target, source in (("Curated without stage", "4 and 6 thermodynamic directions", ["4_thermo", "6_v2_thermo"], cur, draft),
                                                   ("Draft with stage only", "4 directions and GTP/GDP carrier", ["4_thermo", "6_v3_gtp_carrier"], draft, cur),
                                                   ("Draft with stage only", "6 directions and GTP/GDP carrier", ["6_v2_thermo", "6_v3_gtp_carrier"], draft, cur)):
        rids = sorted(set(log[log.stage.isin(stages)].reaction_id))
        m = target.copy()
        n = transfer(m, source, rids)
        rows.append({"variant": variant, "stage": label, "n_reactions": len(rids), **n, **measure(m)})
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "effects.csv", index=False)
    write_provenance(OUT, params={"engine": "cobra static FBA", "draft": DRAFT.name, "curated": CURATED.name, "log": "model/curation/curation_log.csv",
                                  "skipped_actions": sorted(ACTIONS_SKIPPED), "uptake_mmol_g_h": 1.0, "fitted_parameters": "none"})
    print(df[["variant", "stage", "n_reactions", "max_atp_per_glucose", "cycle_tests_failed", "open_exports", "yield_glucose_gg", "yield_fructose_gg"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
