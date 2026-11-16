"""Whether transient glycerol in the batch fermentations follows oxygen limitation. Nothing is fitted.

Timing: glycerol, sugars, online biomass and dissolved O2 at each HPLC sampling time. Online readings are the
ones nearest in time for each replicate.

O2 cap: the Draft and Curated models at the sugar uptake that gives mu_max, with specific O2 uptake capped at a
fraction of the least the model needs at that uptake. Export ranges come from flux variability analysis at the
growth optimum. Sugar uptake is either a bound, as in the dFBA growth problem, or fixed at the kinetic rate, for
a cell that keeps taking up sugar when O2 runs short.

Outputs: results/glycerol_oxygen/{timing.csv, o2_cap.csv, summary.json, provenance.json}

Usage (from analysis/):  python scripts/glycerol_oxygen_test.py
"""

from __future__ import annotations

import json

import cobra
import numpy as np
import pandas as pd

from fbaval import data, simulate
from fbaval.provenance import write_provenance

from fbaval.paths import CURATED, DRAFT, RESULTS

IN = RESULTS / "flask_products"
OUT = RESULTS / "glycerol_oxygen"
MODELS = {"Draft": DRAFT, "Curated": CURATED}
SUGARS = ("glucose", "fructose")
CAPS = (1.0, 0.8, 0.6, 0.4, 0.2)
EXPORTS = {"glycerol": "BOUNDARY_GLe", "ethanol": "BOUNDARY_ETHe", "acetaldehyde": "BOUNDARY_ACALe",
           "acetate": "BOUNDARY_ACe", "pyruvate": "BOUNDARY_PYRe"}
TOL = 1e-6


def timing() -> pd.DataFrame:
    on = pd.read_csv(IN / "online.csv").dropna(subset=["do_pct_as", "biomass_gl"])
    prod = pd.read_csv(IN / "products.csv")
    rows = []
    for cond in data.PRIMARY_VINESTOCK:
        p = prod[prod.condition == cond]
        mean = {c: p[p.compound == c].groupby("time_h").conc_gl.mean() for c in ("glycerol", "glucose", "fructose")}
        peak = mean["glycerol"].max()
        for t in mean["glycerol"].index:
            do, x = [], []
            for _, o in on[on.condition == cond].groupby("flask"):
                i = (o.time_h - t).abs().idxmin()
                do.append(float(o.loc[i, "do_pct_as"]))
                x.append(float(o.loc[i, "biomass_gl"]))
            rows.append({"condition": cond, "time_h": t, "glycerol_gl": mean["glycerol"][t],
                         "share_of_peak": mean["glycerol"][t] / peak if peak > 0 else 0.0,
                         "glucose_gl": mean["glucose"].get(t, 0.0), "fructose_gl": mean["fructose"].get(t, 0.0),
                         "biomass_gl": np.mean(x), "do_mean_pct": np.mean(do), "do_min_pct": min(do), "do_max_pct": max(do)})
    return pd.DataFrame(rows)


def ranges(model, ids: dict, fraction: float = 1.0) -> dict:
    rxns = [model.reactions.get_by_id(r) for r in ids.values() if r in model.reactions]
    fva = cobra.flux_analysis.flux_variability_analysis(model, reaction_list=rxns, fraction_of_optimum=fraction)
    out = {}
    for name, rid in ids.items():
        if rid in fva.index:
            out[f"{name}_min"] = max(float(fva.loc[rid, "minimum"]), 0.0)
            out[f"{name}_max"] = max(float(fva.loc[rid, "maximum"]), 0.0)
    return out


def o2_cap(label: str, model) -> list[dict]:
    q = simulate.vmax_from_mu_max(model, simulate.INPUTS["mu_max_per_h"][0])
    open_exports = None
    if label == "Curated":
        with model:
            simulate.set_sugars(model, {"glucose": q})
            open_exports = {r.id: r.id for r in model.boundary if r.upper_bound > 0 and r.id not in
                            ("BOUNDARY_O2e", "BOUNDARY_CO2e", "BOUNDARY_H2Oe", "BOUNDARY_He", "DRAIN_Biomass")
                            and "C" in next(iter(r.metabolites)).elements}
    rows = []
    for sugar in SUGARS:
        with model:
            simulate.set_sugars(model, {sugar: q})
            model.slim_optimize()
            fva = cobra.flux_analysis.flux_variability_analysis(model, reaction_list=[model.reactions.BOUNDARY_O2e], fraction_of_optimum=1.0)
            need = -float(fva.loc["BOUNDARY_O2e", "maximum"])
        for mode, forced in (("bound", False), ("fixed", True)):
            for cap in CAPS:
                with model:
                    simulate.set_sugars(model, {sugar: q}, o2_cap=need * cap, forced=forced)
                    mu = model.slim_optimize(error_value=float("nan"))
                    row = {"model": label, "sugar": sugar, "uptake": mode, "q_sugar": q, "qo2_need": need, "o2_cap_fraction": cap, "mu": mu}
                    if mu == mu:
                        row.update(ranges(model, EXPORTS))
                        row["glycerol_max_at_99pct_growth"] = ranges(model, {"glycerol": EXPORTS["glycerol"]}, 0.99)["glycerol_max"]
                        if open_exports:
                            r = ranges(model, open_exports)
                            row["forced_carbon_exports"] = ";".join(sorted(k.replace("BOUNDARY_", "") for k in open_exports if r[f"{k}_min"] > TOL))
                rows.append(row)
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    tim = timing()
    tim.to_csv(OUT / "timing.csv", index=False)
    rows = []
    for label, file in MODELS.items():
        rows += o2_cap(label, simulate.load(str(file), normalize=True))
    cap = pd.DataFrame(rows)
    cap.to_csv(OUT / "o2_cap.csv", index=False)

    summary = {"timing": {}, "o2_cap": {}}
    for cond, t in tim.groupby("condition", sort=False):
        if t.glycerol_gl.max() <= 0:
            continue
        pk = t.loc[t.glycerol_gl.idxmax()]
        high = t[(t.do_mean_pct >= 100) & (t.time_h <= pk.time_h)]
        last = high.iloc[-1]
        has = t[t.glucose_gl > 0.05]
        none = t[(t.glucose_gl <= 0.05) & (t.time_h > has.time_h.max())]
        summary["timing"][cond] = {
            "peak_gl": pk.glycerol_gl, "peak_t_h": pk.time_h, "peak_do_mean": pk.do_mean_pct,
            "last_sample_do_ge_100": {"t_h": last.time_h, "glycerol_gl": last.glycerol_gl, "share_of_peak": last.share_of_peak,
                                      "do_mean": last.do_mean_pct, "do_min": last.do_min_pct},
            "biomass_last_sample_with_glucose": has.iloc[-1].biomass_gl,
            "biomass_first_sample_without_glucose": none.iloc[0].biomass_gl if len(none) else None,
            "do_min_any_sample": t.do_min_pct.min()}
    limited = cap[cap.o2_cap_fraction < 1]
    for (label, mode), c in limited.groupby(["model", "uptake"], sort=False):
        summary["o2_cap"][f"{label}, uptake {mode}"] = {
            "glycerol_max_at_optimum": c.glycerol_max.max(), "glycerol_max_at_99pct_growth": c.glycerol_max_at_99pct_growth.max(),
            "ethanol_max": c.ethanol_max.max(), "acetaldehyde_max": c.acetaldehyde_max.max(), "acetate_max": c.acetate_max.max(), "pyruvate_max": c.pyruvate_max.max()}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    write_provenance(OUT, params={"engine": "COBRApy FBA and flux variability analysis", "models": {k: v.name for k, v in MODELS.items()}, "o2_cap_fractions": CAPS,
                                  "mu_max_per_h": simulate.INPUTS["mu_max_per_h"][0], "input": "results/flask_products", "fitted": "nothing"})
    print(tim[tim.condition != "Fru15"].round(3).to_string(index=False))
    show = ["model", "sugar", "uptake", "o2_cap_fraction", "mu", "glycerol_min", "glycerol_max", "glycerol_max_at_99pct_growth", "ethanol_min", "ethanol_max", "acetaldehyde_min", "acetaldehyde_max", "acetate_max", "pyruvate_max"]
    print(cap[show].round(3).to_string(index=False))
    print(cap[cap.model == "Curated"][["sugar", "uptake", "o2_cap_fraction", "forced_carbon_exports"]].to_string(index=False))
    print(json.dumps(summary, indent=2, default=lambda v: round(float(v), 3)))


if __name__ == "__main__":
    main()
