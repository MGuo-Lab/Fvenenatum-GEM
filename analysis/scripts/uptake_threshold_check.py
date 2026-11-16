"""Sensitivity of the batch-fermentation runs to the glucose threshold of sequential uptake, and timing of
infeasible dFBA steps.

1. Threshold. Sequential uptake blocks fructose while glucose is above 0.5 mM. The two mixed conditions are
   rerun for both models with thresholds of 0.05, 0.5 and 5 mM, with x0 refitted per condition at each
   threshold (timing criterion and fit of scripts/flask_dynamics.py). Reported: final biomass, its
   discrepancy, fitted x0 and time-course RMSE.
2. Infeasible steps. For every model and condition in the final setting (sequential uptake, 0.5 mM, x0 fitted
   per condition), the time of the first infeasible step is compared with the time at which all sugars fall
   below 0.01 g/L.

Outputs: results/uptake_threshold_check/{threshold.csv, infeasible.csv, provenance.json}

Usage (from analysis/):  python scripts/uptake_threshold_check.py
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import flask_dynamics as fd
from fbaval import data
from fbaval.paths import RESULTS
from fbaval.provenance import write_provenance

OUT = RESULTS / "uptake_threshold_check"
THRESHOLDS = (0.05, 0.5, 5.0)
MIXED = ("Glc10Fru5", "Glc5Fru10")


def _set_thr(thr: float) -> None:
    fd.THR = thr


def fit_eval(args: tuple[str, str, float]) -> dict:
    key, cond, thr = args
    _set_thr(thr)
    f = fd.fit((key, cond, True))
    row, _ = fd.evaluate((key, "lag_each", cond, f["x0"], True, f["status"]))
    return {**row, "threshold_mM": thr}


def infeasible(args: tuple[str, str]) -> dict:
    key, cond = args
    _set_thr(0.5)
    f = fd.fit((key, cond, True))
    sim = fd.run(key, cond, f["x0"], True)
    sug = sum(sim[c] for c in ("glucose_gl", "fructose_gl") if c in sim)
    t_dep = float(sim.time_h[sug < 0.01].min()) if (sug < 0.01).any() else np.nan
    inf = sim[sim.infeasible.astype(bool)] if "infeasible" in sim else sim.iloc[0:0]
    t_inf = float(inf.time_h.min()) if len(inf) else np.nan
    return {"model": key, "condition": cond, "x0_gl": f["x0"], "t_sugars_depleted_h": t_dep,
            "t_first_infeasible_h": t_inf, "n_infeasible_steps": int(len(inf)), "n_steps": int(len(sim)),
            "x_final_gl": float(sim.biomass_gl.max()), "x_end_gl": float(sim.biomass_gl.iloc[-1])}


def main() -> None:
    with ProcessPoolExecutor(max_workers=12) as ex:
        thr = pd.DataFrame(list(ex.map(fit_eval, [(k, c, t) for k in fd.MODELS for c in MIXED for t in THRESHOLDS])))
        inf = pd.DataFrame(list(ex.map(infeasible, [(k, c) for k in fd.MODELS for c in data.PRIMARY_VINESTOCK])))
    OUT.mkdir(parents=True, exist_ok=True)
    thr.to_csv(OUT / "threshold.csv", index=False)
    inf.to_csv(OUT / "infeasible.csv", index=False)
    write_provenance(OUT, params={"engine": "dfba.run_dfba_multi + cobra FBA (via scripts/flask_dynamics.py)",
                                  "models": {k: v.name for k, v in fd.MODELS.items()}, "dt_h": fd.DT, "t_end": fd.T_END,
                                  "thresholds_mM": THRESHOLDS, "conditions_threshold": MIXED,
                                  "sugar_depleted_below_gl": 0.01,
                                  "fitted_parameters": "x0 per condition at each threshold; timing criterion t50"})
    print(thr[["model", "condition", "threshold_mM", "x0_gl", "x_final_gl", "err_final_pct", "timing_err_h",
               "rmse_x", "rmse_s"]].round(4).to_string(index=False))
    print(inf.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
