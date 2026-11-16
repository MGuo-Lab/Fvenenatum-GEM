"""Batch fermentations: stepwise agreement with biomass and sugar time courses, Draft and Curated.

Settings, each run for both models with nothing else changed:
  base         x0 = 0.01 g/L, simultaneous sugar uptake
  sugar_order  x0 = 0.01 g/L, glucose-first uptake
  lag_only     simultaneous uptake, x0 fitted per condition (no sugar order)
  lag_shared   glucose-first uptake, one x0 fitted across the four conditions
  lag_each     glucose-first uptake, x0 fitted per condition

Lag is represented by the initial biomass x0. It is fitted to one number per condition, the time at
which the measured mean reaches half the mean final dry weight, separately for each model. Literature
kinetics stay fixed: mu_max 0.20 1/h, Km 0.26 mM. Glucose-first uptake blocks other sugars while glucose
is above 0.5 mM. Time-course error is computed over all replicates and all measured time points:
biomass against online biomass, sugars against HPLC.

Outputs: results/flask_dynamics/{metrics.csv, summary.csv, profile.csv, series.csv, provenance.json}

Usage (from analysis/):  python scripts/flask_dynamics.py
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize_scalar

from fbaval import data, simulate
from fbaval.paths import CURATED, DRAFT, RESULTS
from fbaval.provenance import write_provenance

OUT = RESULTS / "flask_dynamics"
MODELS = {"draft": DRAFT, "curated": CURATED}
DT, T_END, X0_BASE = 0.5, 150.0, 0.01
THR = simulate.INPUTS["diauxie_threshold_mM"][0]
GRID = 10.0 ** np.arange(-6.0, -1.49, 0.5)
SETTINGS = ("base", "sugar_order", "lag_only", "lag_shared", "lag_each")
_CACHE: dict[str, object] = {}


def model(key: str):
    if key not in _CACHE:
        _CACHE[key] = simulate.load(str(MODELS[key]), normalize=True)
    return _CACHE[key]


def run(key: str, cond: str, x0: float, glucose_first: bool) -> pd.DataFrame:
    m = model(key)
    vmax = simulate.vmax_from_mu_max(m, simulate.INPUTS["mu_max_per_h"][0])
    km = simulate.INPUTS["km_sugar_mM"][0]
    spec = data.VINESTOCK[cond]
    sugars = {s: spec[f"{s}_gl"] for s in ("glucose", "fructose") if spec[f"{s}_gl"] > 0}
    return simulate.run_dfba(m, sugars, x0, T_END, DT, vmax, km, data.VINESTOCK_OTR_MAX_MMOL_L_H,
                             diauxie_threshold_mM=THR if glucose_first else None)


def t50(sim: pd.DataFrame) -> float:
    return float(np.interp(0.5 * sim.biomass_gl.max(), sim.biomass_gl.cummax(), sim.time_h))


def measured_t50(cond: str) -> float:
    pres = data.load_vinestock_presens(cond).dropna(subset=["biomass_gl"])
    fin = data.load_vinestock_final(cond).biomass_gl
    grid = np.arange(0, pres.time_h.max(), 0.25)
    mean = np.nanmean([np.interp(grid, d.time_h, d.biomass_gl, left=np.nan, right=np.nan)
                       for _, d in pres.groupby("flask")], axis=0)
    return float(grid[np.nanargmax(mean >= 0.5 * fin.mean())])


def fit(args: tuple[str, str, bool]) -> dict:
    """x0 giving the measured t50, from a log grid and a root search; the grid end if none is bracketed."""
    key, cond, gf = args
    t_meas = measured_t50(cond)
    t = np.array([t50(run(key, cond, float(x), gf)) for x in GRID])
    d = t - t_meas
    i = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
    if len(i):
        lx = brentq(lambda v: t50(run(key, cond, 10.0 ** v, gf)) - t_meas,
                    np.log10(GRID[i[0]]), np.log10(GRID[i[0] + 1]), xtol=1e-3)
        status = "fitted"
    else:
        lx = float(np.log10(GRID[int(np.argmin(np.abs(d)))]))
        status = "at grid bound"
    return {"model": key, "condition": cond, "glucose_first": gf, "x0": 10.0 ** lx, "status": status,
            "t_meas": t_meas, "grid_t50": t}


def errors(sim: pd.DataFrame, cond: str) -> dict:
    pres = data.load_vinestock_presens(cond).dropna(subset=["biomass_gl"])
    hplc = data.load_vinestock_hplc(cond)
    ex = np.interp(pres.time_h, sim.time_h, sim.biomass_gl) - pres.biomass_gl.to_numpy()
    es = []
    for s in ("glucose", "fructose"):
        if data.VINESTOCK[cond][f"{s}_gl"] > 0:
            q = hplc[hplc.compound == s]
            es.append(np.interp(q.time_h, sim.time_h, sim[f"{s}_gl"]) - q.conc_gl.to_numpy())
    es = np.concatenate(es)
    return {"rmse_x": float(np.sqrt(np.mean(ex ** 2))), "mae_x": float(np.mean(np.abs(ex))), "n_x": len(ex),
            "rmse_s": float(np.sqrt(np.mean(es ** 2))), "mae_s": float(np.mean(np.abs(es))), "n_s": len(es)}


def evaluate(args: tuple[str, str, str, float, bool, str]) -> tuple[dict, pd.DataFrame]:
    key, setting, cond, x0, gf, status = args
    sim = run(key, cond, x0, gf)
    meas = float(data.load_vinestock_final(cond).biomass_gl.mean())
    x_final = float(sim.biomass_gl.max())
    t_meas = measured_t50(cond)
    row = {"model": key, "setting": setting, "condition": cond, "x0_gl": x0, "status": status, "glucose_first": gf,
           "t50_meas_h": t_meas, "t50_model_h": t50(sim), "timing_err_h": t50(sim) - t_meas, "x_final_gl": x_final,
           "x_final_meas_gl": meas, "err_final_pct": 100 * (x_final / meas - 1), **errors(sim, cond)}
    cols = ["time_h", "biomass_gl"] + [f"{s}_gl" for s in ("glucose", "fructose") if f"{s}_gl" in sim]
    return row, sim[cols].assign(model=key, setting=setting, condition=cond)


def main() -> None:
    conds = list(data.PRIMARY_VINESTOCK)
    jobs = [(k, c, gf) for k in MODELS for c in conds for gf in (True, False)]
    with ProcessPoolExecutor(max_workers=6) as ex:
        fits = list(ex.map(fit, jobs))
    F = {(r["model"], r["condition"], r["glucose_first"]): r for r in fits}

    shared = {}
    for k in MODELS:
        def sse(lx):
            return sum((np.interp(lx, np.log10(GRID), F[(k, c, True)]["grid_t50"]) - F[(k, c, True)]["t_meas"]) ** 2 for c in conds)
        shared[k] = 10.0 ** minimize_scalar(sse, bounds=(-6.0, -1.5), method="bounded").x

    tasks = []
    for k in MODELS:
        for c in conds:
            tasks += [(k, "base", c, X0_BASE, False, "assumed"), (k, "sugar_order", c, X0_BASE, True, "assumed"),
                      (k, "lag_only", c, F[(k, c, False)]["x0"], False, F[(k, c, False)]["status"]),
                      (k, "lag_shared", c, shared[k], True, "fitted, shared"),
                      (k, "lag_each", c, F[(k, c, True)]["x0"], True, F[(k, c, True)]["status"])]
    with ProcessPoolExecutor(max_workers=6) as ex:
        res = list(ex.map(evaluate, tasks))
    metrics = pd.DataFrame([r for r, _ in res])
    series = pd.concat([s for _, s in res], ignore_index=True)
    profile = pd.DataFrame([{"model": r["model"], "condition": r["condition"], "glucose_first": r["glucose_first"], "x0_gl": x, "t50_h": t}
                            for r in fits for x, t in zip(GRID, r["grid_t50"])])
    metrics["abs_err_final_pct"] = metrics.err_final_pct.abs()
    metrics["abs_timing_err_h"] = metrics.timing_err_h.abs()
    keep = ["rmse_x", "mae_x", "rmse_s", "mae_s", "abs_err_final_pct", "abs_timing_err_h"]
    summary = metrics.groupby(["model", "setting"], sort=False)[keep].mean().reset_index()

    OUT.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(OUT / "metrics.csv", index=False)
    summary.to_csv(OUT / "summary.csv", index=False)
    profile.to_csv(OUT / "profile.csv", index=False)
    series.to_csv(OUT / "series.csv", index=False)
    write_provenance(OUT, params={"engine": "dfba.run_dfba_multi + cobra FBA", "models": {k: v.name for k, v in MODELS.items()},
                                  "settings": SETTINGS, "mu_max": simulate.INPUTS["mu_max_per_h"][0],
                                  "km_mM": simulate.INPUTS["km_sugar_mM"][0], "otr_mmol_L_h": data.VINESTOCK_OTR_MAX_MMOL_L_H,
                                  "dt": DT, "t_end": T_END, "x0_base": X0_BASE, "diauxie_threshold_mM": THR,
                                  "fitted_parameters": "x0 only (per condition, or one shared); timing criterion t50"})
    print(summary.round(3).to_string(index=False))
    print(metrics[metrics.setting == "lag_each"][["model", "condition", "x0_gl", "status", "err_final_pct", "rmse_x", "rmse_s"]].round(4).to_string(index=False))
    print(metrics[metrics.setting == "lag_only"][["model", "condition", "x0_gl", "status", "timing_err_h"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
