"""Checks of two analysis choices behind the batch-fermentation comparison, Draft and Curated.

1. dFBA step size. The flask runs use an explicit Euler step of 0.5 h. Each model and condition is rerun
   at 0.5, 0.25 and 0.1 h in two settings: sequential uptake with the assumed x0 of 0.01 g/L (the
   final-biomass test, nothing fitted) and sequential uptake with x0 refitted per condition at each step
   (the final time-course setting). Reported: final biomass, its discrepancy, time to half final biomass,
   fitted x0 and time-course RMSE.
2. Ash. The biomass reaction weighs 1 g/mmol of ash-free dry mass, whereas measured dry weight includes
   ash. With an ash fraction a of measured dry weight, the predicted dry weight is x_pred / (1 - a). The
   discrepancy is recomputed for a = 0.05 (Supplementary Table S2) from the 0.5 h runs.

Fitting, timing criterion and error metrics are those of scripts/flask_dynamics.py.

Outputs: results/methods_checks/{step_size.csv, step_size_summary.csv, ash.csv, ash_summary.csv,
provenance.json}

Usage (from analysis/):  python scripts/methods_checks.py
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import pandas as pd

import flask_dynamics as fd
from fbaval import data
from fbaval.paths import RESULTS
from fbaval.provenance import write_provenance

OUT = RESULTS / "methods_checks"
STEPS = (0.5, 0.25, 0.1)
ASH = 0.05


def _set_dt(dt: float) -> None:
    fd.DT = dt


def fit_at(args: tuple[str, str, float]) -> dict:
    key, cond, dt = args
    _set_dt(dt)
    r = fd.fit((key, cond, True))
    return {"model": key, "condition": cond, "dt_h": dt, "x0_gl": r["x0"], "status": r["status"]}


def evaluate_at(args: tuple[str, str, str, float, float, str]) -> dict:
    key, setting, cond, x0, dt, status = args
    _set_dt(dt)
    row, _ = fd.evaluate((key, setting, cond, x0, True, status))
    return {**row, "dt_h": dt}


def main() -> None:
    conds = list(data.PRIMARY_VINESTOCK)
    with ProcessPoolExecutor(max_workers=24) as ex:
        fits = list(ex.map(fit_at, [(k, c, dt) for k in fd.MODELS for c in conds for dt in STEPS]))
    X0 = {(f["model"], f["condition"], f["dt_h"]): f for f in fits}
    tasks = []
    for k in fd.MODELS:
        for c in conds:
            for dt in STEPS:
                tasks.append((k, "sugar_order", c, fd.X0_BASE, dt, "assumed"))
                f = X0[(k, c, dt)]
                tasks.append((k, "lag_each", c, f["x0_gl"], dt, f["status"]))
    with ProcessPoolExecutor(max_workers=24) as ex:
        rows = list(ex.map(evaluate_at, tasks))
    step = pd.DataFrame(rows)
    step["abs_err_final_pct"] = step.err_final_pct.abs()
    summary = (step.groupby(["model", "setting", "dt_h"], sort=False)
               [["abs_err_final_pct", "rmse_x", "rmse_s", "mae_x", "mae_s"]].mean().reset_index())

    ash = step[(step.setting == "sugar_order") & (step.dt_h == 0.5)][
        ["model", "condition", "x_final_gl", "x_final_meas_gl", "err_final_pct"]].copy()
    ash["ash_fraction"] = ASH
    ash["x_pred_dw_gl"] = ash.x_final_gl / (1 - ASH)
    ash["err_final_ash_pct"] = 100 * (ash.x_pred_dw_gl / ash.x_final_meas_gl - 1)
    ash_summary = ash.assign(abs_err=ash.err_final_pct.abs(), abs_err_ash=ash.err_final_ash_pct.abs()) \
        .groupby("model")[["abs_err", "abs_err_ash"]].mean().reset_index()

    OUT.mkdir(parents=True, exist_ok=True)
    step.to_csv(OUT / "step_size.csv", index=False)
    summary.to_csv(OUT / "step_size_summary.csv", index=False)
    ash.to_csv(OUT / "ash.csv", index=False)
    ash_summary.to_csv(OUT / "ash_summary.csv", index=False)
    write_provenance(OUT, params={"engine": "dfba.run_dfba_multi + cobra FBA (via scripts/flask_dynamics.py)",
                                  "models": {k: v.name for k, v in fd.MODELS.items()}, "dt_h": STEPS, "t_end": fd.T_END, "x0_base": fd.X0_BASE,
                                  "diauxie_threshold_mM": fd.THR, "ash_fraction": ASH,
                                  "fitted_parameters": "x0 per condition at each step size; timing criterion t50"})
    print(summary.round(3).to_string(index=False))
    print(step[["model", "setting", "condition", "dt_h", "x0_gl", "x_final_gl", "err_final_pct", "t50_model_h",
                "rmse_x", "rmse_s"]].round(4).to_string(index=False))
    print(ash.round(3).to_string(index=False))
    print(ash_summary.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
