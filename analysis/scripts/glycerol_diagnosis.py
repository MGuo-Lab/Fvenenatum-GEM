"""Whether transient glycerol in the batch fermentations follows dissolved O2, sugar or biomass (Fig. 4, Results).

Measured data and single FBA solves only.

  1. Clock. The time of the first online reading on the HPLC clock, checked against the temperature dip that follows
     each sampling.
  2. Rates. Glycerol, sugar and growth per interval between HPLC samples, per flask. Rates use the integral of online
     biomass over the interval.
  3. Predictors of glycerol while glucose remains, each fitted on two conditions and tested on the third.
  4. Dissolved O2 against biomass, and the O2 budget: the integral of (200 - DO) dt per flask.
  5. Matched glucose and fructose intervals at the same biomass.
  6. The 16.1 min peak, from results/flask_hplc, against glycerol and against the peak table.
  7. iJL1143: O2 cap with sugar uptake fixed, with electron and carbon exits closed step by step, on both sugars.
     Growth cost of the measured glycerol rate. Cap on cytosolic NADH shuttles.

Steps 1 and 6 read the raw PreSens exports and HPLC data files under FBAVAL_RAW_DIR (see fbaval.paths), which are
not part of this repository. Run scripts/flask_hplc.py first. Outputs: results/glycerol_diagnosis/{clock.csv,
samples.csv, intervals.csv, predictors.csv, do_vs_biomass.csv, o2_budget.csv, matched.csv, peak_16p1.csv,
gem_o2_cap.csv, gem_glycerol_cost.csv, gem_nadh_cap.csv, summary.json, provenance.json}

Usage (from analysis/):  python scripts/glycerol_diagnosis.py
"""

from __future__ import annotations

import json
import logging

import cobra
import numpy as np
import openpyxl
import pandas as pd
from scipy.optimize import least_squares

from fbaval import data, simulate
from fbaval.lcd import chromatogram, read_injection
from fbaval.paths import CURATED, RAW_DIR, RESULTS
from fbaval.provenance import write_provenance

logging.disable(logging.WARNING)
OUT = RESULTS / "glycerol_diagnosis"
RAW = RAW_DIR / "lcd" / "TV" / "Batch Shake Flask Winter 2025-26"
RUNS = {"Glc15": "19012026_Glucose_15gpL", "Glc10Fru5": "02122025GluFru", "Glc5Fru10": "15122025_FruGlu", "Fru15": "15012026_Fru"}
CONDS = data.PRIMARY_VINESTOCK
GLC_CONDS = CONDS[:3]
O2_PER_PCT = data.VINESTOCK_OTR_MAX_MMOL_L_H / 200.0   # mmol/L/h per % air saturation, the paper's kLa 96 /h at 200 %
GLYCEROL_MM = 92.09
HEXOSE_MM = 180.16
GLYCEROL_AREA_PER_GL_UL = 13170.0
Q_GLYCEROL_TEST = 0.55
SKIP = ("BOUNDARY_CO2e", "BOUNDARY_H2Oe", "BOUNDARY_He", "DRAIN_Biomass", "BOUNDARY_O2e")
INORGANIC = ("BOUNDARY_H2SO3e", "BOUNDARY_HNO2", "BOUNDARY_HNO3e", "BOUNDARY_ASe")


def smooth(s: np.ndarray, n: int = 5) -> np.ndarray:
    return pd.Series(s).rolling(n, center=True, min_periods=1).median().to_numpy()


def online(cond: str) -> tuple[pd.DataFrame, float]:
    """Online readings, long form, and the first reading's time after the first HPLC sample."""
    o = data.load_vinestock_presens(cond).rename(columns={"time_h": "t", "do_pct_as": "do", "biomass_gl": "x"})
    return o[["t", "flask", "do", "x"]], float(o.t.min())


def sampling_dips(cond: str) -> pd.DataFrame:
    """Time on the HPLC clock of the temperature dip after each HPLC sample, from the raw PreSens Detail sheets."""
    folder = RAW_DIR / "presens" / data.VINESTOCK[cond]["presens_dir"]
    t0 = pd.to_datetime(pd.read_excel(data.condition_dir(cond) / "Processed_HPLC_Data.xlsx")["Date Time"].iloc[0])
    frames = []
    for f in sorted((folder / "PreSens Data").glob("*.xlsx")):
        ws = openpyxl.load_workbook(f, read_only=True, data_only=True)["Detail"]
        rows = list(ws.iter_rows(values_only=True))
        head = [str(h) for h in rows[0]]
        temp = [i for i, h in enumerate(head) if h.startswith("Temperature")]
        frames.append(pd.DataFrame([(r[1], np.nanmean([float(r[i]) if isinstance(r[i], (int, float)) else np.nan for i in temp])) for r in rows[2:] if r[0] is not None], columns=["abs", "temp"]))
    d = pd.concat(frames).sort_values("abs")
    d["t"] = (pd.to_datetime(d["abs"]) - t0).dt.total_seconds() / 3600.0
    d["drop"] = d.temp.diff()
    return d[d["drop"] < -0.08][["t", "drop"]]


def tables() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    samples, intervals, on = [], [], {}
    for cond in CONDS:
        o, _ = online(cond)
        on[cond] = o
        prod = data.load_vinestock_products(cond)
        hplc = data.load_vinestock_hplc(cond)
        for fl in (1, 2, 3):
            q = prod[(prod.flask == fl) & (prod.compound == "glycerol")].sort_values("time_h")[["time_h", "conc_gl"]].rename(columns={"conc_gl": "glycerol"})
            for sugar in ("glucose", "fructose"):
                s = hplc[(hplc.flask == fl) & (hplc.compound == sugar)].set_index("time_h").conc_gl
                q[sugar] = q.time_h.map(s).fillna(0.0).to_numpy()
            d = o[o.flask == fl].dropna(subset=["do", "x"]).sort_values("t")
            t, do, x = d.t.to_numpy(), smooth(d.do.to_numpy()), smooth(d.x.to_numpy())
            tau = np.concatenate([[0.0], np.cumsum(0.5 * (x[1:] + x[:-1]) * np.diff(t))])
            q["do"], q["x"], q["tau"] = np.interp(q.time_h, t, do), np.interp(q.time_h, t, x), np.interp(q.time_h, t, tau)
            q["glc_ref"] = q.glucose.iloc[1]    # first sample after time zero, because time zero is capped in the workbooks
            q["flask"], q["condition"] = fl, cond
            samples.append(q)
            for (_, a), (_, b) in zip(q.iloc[:-1].iterrows(), q.iloc[1:].iterrows()):
                m = (t >= a.time_h) & (t <= b.time_h)
                if m.sum() < 2:
                    continue
                ix = np.trapezoid(x[m], t[m])
                intervals.append({"condition": cond, "flask": fl, "t_start": a.time_h, "t_end": b.time_h, "x_start": a.x, "x_end": b.x,
                                  "do_mean": do[m].mean(), "do_min": do[m].min(), "glycerol_start": a.glycerol, "glycerol_end": b.glycerol,
                                  "glucose_start": a.glucose, "glucose_end": b.glucose, "fructose_start": a.fructose, "fructose_end": b.fructose,
                                  "q_glycerol": (b.glycerol - a.glycerol) / GLYCEROL_MM * 1000 / ix if ix > 0 else np.nan,
                                  "q_glucose": (a.glucose - b.glucose) / HEXOSE_MM * 1000 / ix if ix > 0 else np.nan,
                                  "q_fructose": (a.fructose - b.fructose) / HEXOSE_MM * 1000 / ix if ix > 0 else np.nan,
                                  "mu": np.log(b.x / a.x) / (b.time_h - a.time_h) if a.x > 0 and b.x > 0 else np.nan})
    return pd.concat(samples, ignore_index=True), pd.DataFrame(intervals), on


def saturating(p, z):
    return p[0] * (1 - np.exp(-p[1] * z))


def predictors(samples: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Glucose-phase samples. Each covariate is fitted on two conditions and tested on the third."""
    a = samples[samples.condition.isin(GLC_CONDS) & (samples.time_h >= 20) & (samples.glucose > 0)].copy()
    a["do_drop_200"] = (200 - a.do).clip(lower=0)
    a["glucose_used"] = (a.glc_ref - a.glucose).clip(lower=0)
    covariates = {"biomass": "x", "biomass-time": "tau", "glucose used": "glucose_used", "O2 drop below 200 %": "do_drop_200"}
    rows = []
    for name, col in covariates.items():
        for sens, sub in (("all flasks", a), ("without Glc 10 + Fru 5 flask 1", a[~((a.condition == "Glc10Fru5") & (a.flask == 1))])):
            for held in GLC_CONDS:
                tr, te = sub[sub.condition != held], sub[sub.condition == held]
                fit = least_squares(lambda p: saturating(p, tr[col].to_numpy()) - tr.glycerol.to_numpy(), [0.7, 0.1], bounds=(0, np.inf))
                pred = saturating(fit.x, te[col].to_numpy())
                late = te.time_h >= te.time_h.max() - 1e-6
                rows.append({"covariate": name, "flasks": sens, "held_out": held, "plateau_gl": fit.x[0], "rate": fit.x[1],
                             "rmse_gl": float(np.sqrt(np.mean((pred - te.glycerol) ** 2))), "bias_gl": float(np.mean(pred - te.glycerol)),
                             "bias_last_sample_gl": float(np.mean(pred[late] - te.glycerol[late]))})
    return pd.DataFrame(rows), a


def do_vs_biomass(on: dict) -> pd.DataFrame:
    rows = []
    for cond, o in on.items():
        for fl, d in o.groupby("flask"):
            d = d.dropna(subset=["do", "x"]).sort_values("t")
            x, do = smooth(d.x.to_numpy(), 7), smooth(d.do.to_numpy())
            for target in (0.5, 1, 1.5, 2, 3, 4, 5, 6):
                i = np.where(x >= target)[0]
                rows.append({"condition": cond, "flask": fl, "x_gl": target, "do": do[i[0]] if len(i) else np.nan, "t_h": d.t.to_numpy()[i[0]] if len(i) else np.nan})
    return pd.DataFrame(rows)


def o2_budget(on: dict) -> pd.DataFrame:
    rows = []
    for cond, o in on.items():
        for fl, d in o.groupby("flask"):
            d = d.dropna(subset=["do"]).sort_values("t")
            d = d[d.t >= 24]
            t, do = d.t.to_numpy(), d.do.to_numpy()
            pre = float(np.median(do[t <= 26.5]))
            end = float(np.median(do[t >= t.max() - 1]))
            base = pre + (end - pre) * (t - t[0]) / (t[-1] - t[0])
            rows.append({"condition": cond, "flask": fl, "do_min": do.min(), "do_end": end,
                         "o2_mmol_l_cstar_200": np.trapezoid(np.maximum(200 - do, 0), t) * O2_PER_PCT,
                         "o2_mmol_l_cstar_linear": np.trapezoid(np.maximum(base - do, 0), t) * O2_PER_PCT})
    return pd.DataFrame(rows)


def matched(intervals: pd.DataFrame) -> pd.DataFrame:
    m = intervals[(intervals.x_start >= 0.4) & (intervals.x_end <= 1.9) & (intervals.t_start > 25) & intervals.condition.isin(["Glc15", "Fru15"])]
    return m.groupby("condition").agg(n=("mu", "size"), mu_mean=("mu", "mean"), mu_sd=("mu", "std"), do_mean=("do_mean", "mean"),
                                      q_glycerol_mean=("q_glycerol", "mean"), q_glycerol_sd=("q_glycerol", "std")).reset_index()


def poly_area(t, y, lo, hi, wlo, whi, deg) -> float:
    side = ((t >= wlo) & (t < lo)) | ((t > hi) & (t <= whi))
    mid = (t >= lo) & (t <= hi)
    c = np.polyfit(t[side] - lo, y[side], deg)
    return float(np.trapezoid(y[mid] - np.polyval(c, t[mid] - lo), t[mid] * 60))


def peak_16p1() -> pd.DataFrame:
    """The 16.1 min area of results/flask_hplc with glycerol on the same scale, the table area per uL, and glycerol
    read above a fitted baseline, which also sees the signal the peak table leaves out."""
    unid = pd.read_csv(RESULTS / "flask_hplc" / "unidentified.csv")
    rows = []
    for cond, folder in RUNS.items():
        prod = data.load_vinestock_products(cond).query("compound == 'glycerol'").set_index(["time_h", "flask"]).area
        for _, u in unid[unid.condition == cond].iterrows():
            f = next((RAW / folder).rglob(u.file))
            inj = read_injection(f)
            t, y = (np.asarray(v, float) for v in chromatogram(f, inj.rid))
            rows.append({"condition": cond, "tp": u.tp, "flask": u.flask, "time_h": u.time_h, "file": u.file, "area_16p1": u.peak_16p1_area_per_ul,
                         "table_area_per_ul": u.peak_16p1_table_area_per_ul, "glycerol_area_per_ul": prod.get((u.time_h, u.flask), 0.0) / u.inj_uL,
                         "glycerol_window": poly_area(t, y, 19.45, 20.25, 18.9, 20.8, 2) / inj.inj_vol_uL})
    d = pd.DataFrame(rows)
    d["ratio_to_glycerol"] = d.area_16p1 / d.glycerol_area_per_ul.where(d.glycerol_area_per_ul > 1000)
    return d


def gem() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    m = simulate.load(str(CURATED), normalize=True)
    carbon = lambda r: "C" in next(iter(r.metabolites)).elements
    def exits(s):
        return {r.id.replace("BOUNDARY_", ""): round(float(s.fluxes[r.id]), 3) for r in m.boundary if s.fluxes[r.id] > 1e-4 and r.id not in SKIP}
    caps, cost, nadh = [], [], []
    modes = {"A all exits": (False, None), "B sulfite and nitrite closed": (True, None),
             "C B, glycerol the only carbon exit": (True, ("BOUNDARY_GLe",)), "D B, glycerol or mannitol": (True, ("BOUNDARY_GLe", "BOUNDARY_MNTe"))}
    for mu0 in (0.20, 0.27):
        q = simulate.vmax_from_mu_max(m, mu0)
        for sugar in ("glucose", "fructose"):
            with m:
                simulate.set_sugars(m, {sugar: q})
                m.slim_optimize()
                need = -float(cobra.flux_analysis.flux_variability_analysis(m, reaction_list=[m.reactions.BOUNDARY_O2e], fraction_of_optimum=1.0).loc["BOUNDARY_O2e", "maximum"])
            for cap in (1.0, 0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.4, 0.2):
                for mode, (close_inorganic, keep) in modes.items():
                    with m:
                        simulate.set_sugars(m, {sugar: q}, o2_cap=need * cap, forced=True)
                        if close_inorganic:
                            for i in INORGANIC:
                                m.reactions.get_by_id(i).upper_bound = 0
                        if keep:
                            for r in m.boundary:
                                if r.upper_bound > 0 and r.id not in SKIP and r.id not in keep and carbon(r):
                                    r.upper_bound = 0
                        s = m.optimize()
                        ok = s.status == "optimal"
                        caps.append({"mu_uptake": mu0, "sugar": sugar, "q_sugar": q, "o2_cap": cap, "mode": mode, "mu": s.objective_value if ok else np.nan,
                                     "glycerol": float(s.fluxes["BOUNDARY_GLe"]) if ok else np.nan, "exits": json.dumps(exits(s)) if ok else "infeasible"})
            with m:
                simulate.set_sugars(m, {sugar: q})
                base = m.slim_optimize()
                m.reactions.BOUNDARY_GLe.lower_bound = Q_GLYCEROL_TEST
                cost.append({"mu_uptake": mu0, "sugar": sugar, "q_glycerol": Q_GLYCEROL_TEST, "mu_free": base, "mu_forced": m.slim_optimize()})
            for cap in (None, 3, 2, 1, 0.5, 0):
                with m:
                    simulate.set_sugars(m, {sugar: q}, forced=True)
                    m.reactions.r262.upper_bound = 0     # alanine dehydrogenase as an NADH sink off
                    if cap is not None:
                        m.reactions.r45a.lower_bound = -cap   # cytosolic malate dehydrogenase, NADH use capped
                    s = m.optimize()
                    nadh.append({"mu_uptake": mu0, "sugar": sugar, "shuttle_cap": cap, "mu": s.objective_value, "glycerol": float(s.fluxes["BOUNDARY_GLe"]),
                                 "sorbitol_route_r197": float(s.fluxes["r197"]), "exits": json.dumps(exits(s))})
    return pd.DataFrame(caps), pd.DataFrame(cost), pd.DataFrame(nadh)


def caps_need_ratio(caps: pd.DataFrame) -> float:
    """O2 per sugar in the model on glucose without a cap, mol/mol, from the fluxes at mu 0.20 (qO2 3.99 at q 2.17)."""
    m = simulate.load(str(CURATED), normalize=True)
    q = float(caps[caps.mu_uptake == 0.20].q_sugar.iloc[0])
    with m:
        simulate.set_sugars(m, {"glucose": q})
        return -float(m.optimize().fluxes["BOUNDARY_O2e"]) / q


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    samples, intervals, on = tables()
    pred, _ = predictors(samples)
    dox = do_vs_biomass(on)
    budget = o2_budget(on)
    match = matched(intervals)
    p16 = peak_16p1()
    caps, cost, nadh = gem()

    clock = []
    for c in CONDS:
        first = online(c)[1]
        taps = data.load_vinestock_products(c).time_h.drop_duplicates().sort_values().to_numpy()
        dips = sampling_dips(c)
        lag = [float(dips.t[dips.t >= t - 0.02].min() - t) for t in taps[1:] if (dips.t >= t - 0.02).any()]
        clock.append({"condition": c, "first_online_reading_h": first, "median_dip_after_sample_h": float(np.median(lag)), "max_dip_after_sample_h": float(np.max(lag))})
    clock = pd.DataFrame(clock)

    for name, table in {"clock": clock, "samples": samples, "intervals": intervals, "predictors": pred, "do_vs_biomass": dox, "o2_budget": budget, "matched": match,
                        "peak_16p1": p16, "gem_o2_cap": caps, "gem_glycerol_cost": cost, "gem_nadh_cap": nadh}.items():
        table.to_csv(OUT / f"{name}.csv", index=False)

    o2_need = 1000 * 15.0 / HEXOSE_MM * caps_need_ratio(caps)

    steady = intervals[intervals.condition.isin(GLC_CONDS) & (intervals.glucose_end > 2) & (intervals.x_start > 0.15) & (intervals.t_start > 25)]
    c_first = intervals[(intervals.condition == "Glc5Fru10") & (intervals.t_start > 20) & (intervals.t_end < 38)]
    late = intervals[(intervals.condition == "Glc15") & (intervals.t_start == 46)]
    after = intervals[(intervals.glucose_start == 0) & (intervals.glycerol_start > 0.05) & (intervals.fructose_start > 0)]
    fru = samples[samples.condition == "Fru15"]
    with_all = pred[pred.flasks == "all flasks"]
    # the window above a fitted baseline reads lower than the peak table, so the fructose bound is put on the table's scale
    known = p16[(p16.condition == "Glc15") & (p16.glycerol_area_per_ul > 1000)]
    window_to_table = float((known.glycerol_area_per_ul / known.glycerol_window).median())
    summary = {
        "clock_first_online_reading_h": dict(zip(clock.condition, clock.first_online_reading_h.round(3))),
        "q_glycerol_glucose_above_2gl": {"n": int(len(steady)), "median": steady.q_glycerol.median(), "q25": steady.q_glycerol.quantile(0.25), "q75": steady.q_glycerol.quantile(0.75),
                                         "r_with_do": float(np.corrcoef(steady.q_glycerol, steady.do_mean)[0, 1]), "do_range": [steady.do_mean.min(), steady.do_mean.max()]},
        "q_glycerol_Glc5Fru10_first_interval": c_first.q_glycerol.round(3).tolist(),
        "q_glycerol_Glc15_46_to_48h": {"q": late.q_glycerol.round(3).tolist(), "do_mean": late.do_mean.round(1).tolist()},
        "q_glycerol_after_glucose_gone": after[["condition", "flask", "t_start", "t_end", "q_glycerol", "q_fructose", "do_mean"]].round(3).to_dict("records"),
        "matched_x_0p4_to_1p9": match.round(3).to_dict("records"),
        "predictor_bias_held_out_Glc5Fru10": with_all[with_all.held_out == "Glc5Fru10"].set_index("covariate").bias_gl.round(3).to_dict(),
        "predictor_mean_rmse": with_all.groupby("covariate").rmse_gl.mean().round(3).to_dict(),
        "predictor_mean_rmse_without_flask1": pred[pred.flasks != "all flasks"].groupby("covariate").rmse_gl.mean().round(3).to_dict(),
        "o2_budget_cstar200": budget.groupby("condition", sort=False).o2_mmol_l_cstar_200.apply(lambda s: s.round(0).tolist()).to_dict(),
        "o2_budget_cstar_linear": budget.groupby("condition", sort=False).o2_mmol_l_cstar_linear.apply(lambda s: s.round(0).tolist()).to_dict(),
        "window_to_table_scale_glycerol": window_to_table,
        "fructose_glycerol_window_upper_bound_gl": float(p16[p16.condition == "Fru15"].glycerol_window.max() * window_to_table / GLYCEROL_AREA_PER_GL_UL),
        "fru15_do_range_while_growing": [float(fru[(fru.time_h >= 45) & (fru.time_h <= 65)].do.min()), float(fru[(fru.time_h >= 45) & (fru.time_h <= 65)].do.max())],
        "ratio_16p1_to_glycerol_Glc15": {"mean": p16[p16.condition == "Glc15"].ratio_to_glycerol.mean(), "sd": p16[p16.condition == "Glc15"].ratio_to_glycerol.std(), "n": int(p16[p16.condition == "Glc15"].ratio_to_glycerol.notna().sum())},
        "gem_q_sugar": caps.groupby("mu_uptake").q_sugar.first().round(3).to_dict(),
        "o2_required_model_mmol_l_15gl": float(o2_need),
        "gem_glycerol_at_cap_0p9_mode_C": caps[(caps.o2_cap == 0.9) & (caps["mode"].str.startswith("C"))][["mu_uptake", "sugar", "glycerol"]].round(3).to_dict("records"),
        "gem_growth_loss_pct_at_q_0p55": {f"{r.sugar} {r.mu_uptake}": round(100 * (1 - r.mu_forced / r.mu_free), 1) for r in cost.itertuples()},
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=float))
    write_provenance(OUT, params={"engine": "COBRApy FBA and flux variability analysis (iJL1143), and measured data", "model": CURATED.name, "fitted": "saturating curves in predictors.csv only; nothing else",
                                  "kla_note": "O2 scale is the paper's 46.9 mmol/L/h at 200 %", "inputs": ["results/flask_hplc/unidentified.csv", "raw PreSens files", "raw HPLC .lcd files"]})
    pd.set_option("display.width", 200, "display.max_columns", 40, "display.max_rows", 200)
    print(clock.round(3).to_string(index=False))
    print(json.dumps(summary, indent=1, default=lambda v: round(float(v), 3))[:7000])
    print(with_all.pivot_table(index="covariate", columns="held_out", values=["rmse_gl", "bias_gl"]).round(3).to_string())
    print(pred[pred.flasks != "all flasks"].pivot_table(index="covariate", columns="held_out", values="bias_gl").round(3).to_string())
    print(p16[p16.glycerol_area_per_ul > 1000][["condition", "tp", "flask", "table_area_per_ul", "area_16p1", "glycerol_area_per_ul", "ratio_to_glycerol"]].round(3).to_string(index=False))
    print(dox.pivot_table(index=["condition", "flask"], columns="x_gl", values="do", sort=False).round(0).to_string())


if __name__ == "__main__":
    main()
