"""Batch fermentations: initial biomass x0 as the lag parameter of the Curated dFBA (Supplementary Note 1).

The dFBA has no lag phase; the a priori runs start from an assumed x0 of 0.01 g/L. Here x0 is fitted
instead, per flask condition and as one value shared by all four, so that the Curated run reaches half
its final biomass when the measured mean does. A lower x0 only delays growth: final biomass, the
validation metric, hardly depends on it. Nothing else is fitted.

The fits use glucose-first uptake in the mixed flasks (every other sugar blocked while glucose is
above the 0.5 mM threshold). With simultaneous uptake, as in the a priori runs, each sugar has its own full uptake capacity, the mixed flasks grow at about twice
mu_max and no x0 on the grid reaches the measured timing; that profile is kept for comparison
(rule "simultaneous" in profile.csv).

Inputs per flask: measured sugars; literature mu_max 0.20 1/h, Km 0.26 mM; measured OTR limit
46.9 mmol/L/h; dt 0.5 h, 150 h horizon. Measured timing: mean PreSens biomass of the three flasks
reaching half the mean final dry weight.

Outputs: results/lag_x0_fit/{profile.csv, fits.csv, series.csv, provenance.json}
  profile.csv  t50 of the Curated run on a log grid of x0, per condition and uptake rule
  fits.csv     measured t50, fitted x0 per condition, shared x0 and the timing error at each
  series.csv   Curated time courses: a priori (x0 = 0.01, simultaneous uptake), and
               glucose-first at the fitted x0 and at the shared x0

Usage (from analysis/):  python scripts/lag_x0_fit.py
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize_scalar

from fbaval import data, simulate
from fbaval.paths import CURATED, RESULTS
from fbaval.provenance import write_provenance

MODEL = CURATED
OUT = RESULTS / "lag_x0_fit"
DT, T_END, X0_APRIORI = 0.5, 150.0, 0.01
THR = simulate.INPUTS["diauxie_threshold_mM"][0]      # glucose-first uptake
GRID = 10.0 ** np.arange(-6.0, -1.49, 0.5)            # 1e-6 ... 3.2e-2 g/L
_M = None


def model():
    global _M
    if _M is None:
        _M = simulate.load(str(MODEL), normalize=True)
    return _M


def measured_t50(cond: str) -> float:
    """Time at which the mean PreSens biomass of the three flasks reaches half the mean final dry weight."""
    pres = data.load_vinestock_presens(cond).dropna(subset=["biomass_gl"])
    fin = data.load_vinestock_final(cond).biomass_gl
    grid = np.arange(0, pres.time_h.max(), 0.25)
    mean = np.nanmean([np.interp(grid, d.time_h, d.biomass_gl, left=np.nan, right=np.nan)
                       for _, d in pres.groupby("flask")], axis=0)
    return float(grid[np.nanargmax(mean >= 0.5 * fin.mean())])


def run(cond: str, x0: float, glucose_first: bool = True) -> pd.DataFrame:
    m = model()
    vmax = simulate.vmax_from_mu_max(m, simulate.INPUTS["mu_max_per_h"][0])
    km = simulate.INPUTS["km_sugar_mM"][0]
    spec = data.VINESTOCK[cond]
    sugars = {s: spec[f"{s}_gl"] for s in ("glucose", "fructose") if spec[f"{s}_gl"] > 0}
    return simulate.run_dfba(m, sugars, x0, T_END, DT, vmax, km, data.VINESTOCK_OTR_MAX_MMOL_L_H,
                             diauxie_threshold_mM=THR if glucose_first else None)


def mixed(cond: str) -> bool:
    spec = data.VINESTOCK[cond]
    return spec["glucose_gl"] > 0 and spec["fructose_gl"] > 0


def t50(sim: pd.DataFrame) -> float:
    """Time at which the run reaches half its own final biomass (as in the Fig. 3 shift)."""
    return float(np.interp(0.5 * sim.biomass_gl.max(), sim.biomass_gl.cummax(), sim.time_h))


def fit_condition(cond: str) -> dict:
    t_meas = measured_t50(cond)
    prof = []
    for rule in ("glucose_first", "simultaneous") if mixed(cond) else ("glucose_first",):
        for x0 in GRID:
            sim = run(cond, float(x0), rule == "glucose_first")
            prof.append({"condition": cond, "rule": rule, "x0_gl": float(x0), "t50_h": t50(sim),
                         "x_final_gl": float(sim.biomass_gl.max())})
    prof = pd.DataFrame(prof)
    pf = prof[prof.rule == "glucose_first"].reset_index(drop=True)
    # t50 falls monotonically with x0; solve t50(x0) = t_meas on log10(x0) inside the bracketing grid interval
    f = lambda lx: t50(run(cond, 10.0 ** lx)) - t_meas
    d = pf.t50_h.to_numpy() - t_meas
    i = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
    if len(i):
        lx = brentq(f, np.log10(pf.x0_gl[i[0]]), np.log10(pf.x0_gl[i[0] + 1]), xtol=1e-3)
        status = "fitted"
    else:                                                 # outside the grid: take the closer end
        lx = float(np.log10(pf.x0_gl[int(np.argmin(np.abs(d)))]))
        status = "at grid bound"
    sim_fit = run(cond, 10.0 ** lx)
    sim_apr = run(cond, X0_APRIORI, glucose_first=False)          # the a priori run of Fig. 3
    return {"condition": cond, "t50_meas_h": t_meas, "x0_fit_gl": 10.0 ** lx, "status": status,
            "t50_fit_h": t50(sim_fit), "t50_apriori_h": t50(sim_apr),
            "x_final_fit_gl": float(sim_fit.biomass_gl.max()), "x_final_apriori_gl": float(sim_apr.biomass_gl.max()),
            "x_final_meas_gl": float(data.load_vinestock_final(cond).biomass_gl.mean()),
            "profile": prof, "sim_fit": sim_fit, "sim_apriori": sim_apr}


def shared_run(args) -> tuple[str, pd.DataFrame]:
    cond, x0 = args
    return cond, run(cond, x0)


def main() -> None:
    conds = list(data.PRIMARY_VINESTOCK)
    with ProcessPoolExecutor(max_workers=4) as ex:
        res = list(ex.map(fit_condition, conds))
    prof = pd.concat([r["profile"] for r in res], ignore_index=True)

    # one x0 shared by all four flasks: least squares on the timing error, by interpolating each
    # condition's profile in log10(x0) (t50 is close to linear in log x0), then confirmed by full runs
    gf = prof[prof.rule == "glucose_first"]

    def sse(lx):
        return sum((np.interp(lx, np.log10(g.x0_gl), g.t50_h) - r["t50_meas_h"]) ** 2
                   for r, (_, g) in zip(res, gf.groupby("condition", sort=False)))
    lx_shared = minimize_scalar(sse, bounds=(-6.0, -1.5), method="bounded").x
    with ProcessPoolExecutor(max_workers=4) as ex:
        shared = dict(ex.map(shared_run, [(c, 10.0 ** lx_shared) for c in conds]))

    rows, series = [], []
    for r in res:
        c = r["condition"]
        s_sh = shared[c]
        rows.append({k: r[k] for k in ("condition", "t50_meas_h", "x0_fit_gl", "status", "t50_fit_h", "t50_apriori_h",
                                       "x_final_fit_gl", "x_final_apriori_gl", "x_final_meas_gl")}
                    | {"x0_shared_gl": 10.0 ** lx_shared, "t50_shared_h": t50(s_sh),
                       "timing_err_apriori_h": r["t50_apriori_h"] - r["t50_meas_h"],
                       "timing_err_fit_h": r["t50_fit_h"] - r["t50_meas_h"],
                       "timing_err_shared_h": t50(s_sh) - r["t50_meas_h"],
                       "x_final_shared_gl": float(s_sh.biomass_gl.max())})
        cols = ["time_h", "biomass_gl"] + [f"{s}_gl" for s in ("glucose", "fructose") if f"{s}_gl" in s_sh]
        for run_lab, sim in (("apriori", r["sim_apriori"]), ("fit", r["sim_fit"]), ("shared", s_sh)):
            series.append(sim[cols].assign(condition=c, run=run_lab))
    fits = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    prof.to_csv(OUT / "profile.csv", index=False)
    fits.to_csv(OUT / "fits.csv", index=False)
    pd.concat(series, ignore_index=True).to_csv(OUT / "series.csv", index=False)
    write_provenance(OUT, params={"engine": "dfba.run_dfba_multi + cobra FBA",
                                  "model": str(MODEL.name), "mu_max": simulate.INPUTS["mu_max_per_h"][0],
                                  "km_mM": simulate.INPUTS["km_sugar_mM"][0],
                                  "otr_mmol_L_h": data.VINESTOCK_OTR_MAX_MMOL_L_H, "dt": DT, "t_end": T_END,
                                  "x0_apriori": X0_APRIORI, "x0_grid": [float(x) for x in GRID],
                                  "uptake_rule": f"fits: glucose-first in mixed flasks (threshold {THR} mM); "
                                                 "a priori: simultaneous",
                                  "fitted_parameters": "x0 only (per condition, and one shared); timing criterion t50",
                                  "conditions": conds})
    pd.set_option("display.width", 200)
    print(fits.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
