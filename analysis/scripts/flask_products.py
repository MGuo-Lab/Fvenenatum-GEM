"""Batch fermentations: measured products, sugars, dissolved O2 and yields per condition. No model is run.

Glycerol, ethanol and acetate come from peak areas and each compound's own calibration divisor
(fbaval.data.load_vinestock_products). The workbook "lactate" is left out: its area is an unidentified peak
at 16.1 min, reported as peak area per uL from results/flask_hplc. Run scripts/flask_hplc.py first. Yield follows the source paper: final dry mass / (initial volume
minus all sampled volumes) / nominal initial sugar. The time of sugar depletion extrapolates the last
two HPLC points above 0.5 g/L to zero, limited by the next sample.

Outputs: results/flask_products/{products.csv, online.csv, replicates.csv, summary.csv, provenance.json}

Usage (from analysis/):  python scripts/flask_products.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from fbaval import data
from fbaval.paths import RESULTS
from fbaval.provenance import write_provenance

OUT = RESULTS / "flask_products"
NOMINAL_SUGAR_GL = 15.0


def depletion_time(s: pd.DataFrame) -> float:
    above = s[s.conc_gl > 0.5].sort_values("time_h")
    last, t_last = above.conc_gl.iloc[-1], above.time_h.iloc[-1]
    two = above.tail(2)
    slope = (two.conc_gl.iloc[1] - two.conc_gl.iloc[0]) / (two.time_h.iloc[1] - two.time_h.iloc[0])
    t_next = s[s.time_h > t_last].time_h.min()
    return float(t_next if slope >= 0 else min(max(t_last + last / -slope, t_last), t_next))


def main() -> None:
    products, online, reps = [], [], []
    unid = pd.read_csv(RESULTS / "flask_hplc" / "unidentified.csv")
    for c in data.PRIMARY_VINESTOCK:
        prod = data.load_vinestock_products(c).assign(condition=c)
        prod = prod[prod.compound != "lactate"]
        hplc = data.load_vinestock_hplc(c)
        sug = hplc[hplc.compound.isin(["glucose", "fructose"])]
        total = sug.groupby(["flask", "time_h"]).conc_gl.sum().reset_index()
        pres = data.load_vinestock_presens(c)
        fin = data.load_vinestock_final(c)
        products.append(prod)
        products.append(sug.assign(condition=c, area=np.nan, conc_stored_gl=sug.conc_gl))
        online.append(pres[["time_h", "flask", "biomass_gl", "do_pct_as"]].assign(condition=c))
        for fl in (1, 2, 3):
            s = total[total.flask == fl].sort_values("time_h")
            d = pres[pres.flask == fl].sort_values("time_h")
            b = d.dropna(subset=["biomass_gl"])
            do = d.dropna(subset=["do_pct_as"])
            t_dep = depletion_time(s)
            row = {"condition": c, "flask": fl, "x_final_gl": float(fin.biomass_gl[fl - 1]),
                   "yield_gg": float(fin.biomass_gl[fl - 1]) / NOMINAL_SUGAR_GL, "t_depletion_h": t_dep,
                   "t_end_h": float(b.time_h.max()), "x_online_at_depletion_gl": float(np.interp(t_dep, b.time_h, b.biomass_gl.rolling(5, center=True, min_periods=1).median())),
                   "do_min_pct": float(do.do_pct_as.min())}
            if data.VINESTOCK[c]["glucose_gl"] > 0:
                g_s = sug[(sug.flask == fl) & (sug.compound == "glucose")].sort_values("time_h")
                row["t_glucose_depletion_h"] = depletion_time(g_s)
                row["t_last_with_glucose_h"] = float(g_s[g_s.conc_gl > 0.05].time_h.max())
                row["t_first_without_glucose_h"] = float(g_s[g_s.time_h > row["t_last_with_glucose_h"]].time_h.min())
                row["x_online_at_glucose_depletion_gl"] = float(np.interp(row["t_glucose_depletion_h"], b.time_h, b.biomass_gl.rolling(5, center=True, min_periods=1).median()))
                row["do_at_glucose_depletion_pct"] = float(np.interp(row["t_glucose_depletion_h"], do.time_h, do.do_pct_as.rolling(5, center=True, min_periods=1).median()))
            u = unid[(unid.condition == c) & (unid.flask == fl)].sort_values("time_h")
            row["peak_16p1_max_area_per_ul"] = float(u.peak_16p1_area_per_ul.max())
            row["peak_16p1_max_t_h"] = float(u.time_h.iloc[int(u.peak_16p1_area_per_ul.argmax())])
            row["peak_16p1_last_area_per_ul"] = float(u.peak_16p1_area_per_ul.iloc[-1])
            for comp in ("glycerol", "ethanol", "acetate"):
                q = prod[(prod.compound == comp) & (prod.flask == fl)].sort_values("time_h")
                row[f"{comp}_peak_gl"] = float(q.conc_gl.max())
                row[f"{comp}_peak_t_h"] = float(q.time_h.iloc[int(q.conc_gl.argmax())])
                row[f"{comp}_last_gl"] = float(q.conc_gl.iloc[-1])
                row[f"{comp}_stored_peak_gl"] = float(q.conc_stored_gl.max())
            g = prod[(prod.compound == "glycerol") & (prod.flask == fl)].sort_values("time_h")
            first = g[g.conc_gl > 0.05].time_h.min()
            row["do_while_glycerol_rises_pct"] = float(np.nanmin(np.interp(g[(g.time_h >= first) & (g.time_h <= row["glycerol_peak_t_h"])].time_h, do.time_h, do.do_pct_as))) if np.isfinite(first) else np.nan
            reps.append(row)
    reps = pd.DataFrame(reps)
    numeric = [k for k in reps if k not in ("condition", "flask")]
    summary = reps.groupby("condition", sort=False)[numeric].agg(["mean", "std", "min", "max"])
    summary.columns = [f"{a}_{b}" for a, b in summary.columns]
    OUT.mkdir(parents=True, exist_ok=True)
    pd.concat(products, ignore_index=True).to_csv(OUT / "products.csv", index=False)
    pd.concat(online, ignore_index=True).to_csv(OUT / "online.csv", index=False)
    reps.to_csv(OUT / "replicates.csv", index=False)
    summary.reset_index().to_csv(OUT / "summary.csv", index=False)
    write_provenance(OUT, params={"engine": "none (measured data only)", "nominal_sugar_gl": NOMINAL_SUGAR_GL,
                                  "workbook_lactate": "left out; unidentified 16.1 min peak, see results/flask_hplc",
                                  "fitted_parameters": "none"})
    print(reps.dropna(subset=["t_glucose_depletion_h"])[["condition", "flask", "t_last_with_glucose_h", "t_first_without_glucose_h", "t_glucose_depletion_h", "x_online_at_glucose_depletion_gl", "do_at_glucose_depletion_pct", "t_depletion_h"]].round(2).to_string(index=False))
    show = ["yield_gg_mean", "yield_gg_std", "glycerol_peak_gl_mean", "peak_16p1_max_area_per_ul_mean", "peak_16p1_last_area_per_ul_max", "ethanol_peak_gl_max", "acetate_peak_gl_max",
            "glycerol_last_gl_max", "do_min_pct_min", "do_min_pct_max", "do_while_glycerol_rises_pct_mean", "t_depletion_h_mean", "t_end_h_mean", "x_online_at_depletion_gl_mean"]
    print(summary[show].round(3).T.to_string())


if __name__ == "__main__":
    main()
