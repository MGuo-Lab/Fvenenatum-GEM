"""How much the assumed sugar Km matters in the batch fermentations.

The model uses one Km for both sugars (0.26 mM, assumed). The A. niger high-affinity transporter MstA has
Km 0.025 mM for glucose and 4.5 mM for fructose (van Kuyk et al. 2004), so these two values bracket it.
Initial biomass stays at each model's per-condition fit, glucose-first uptake stays on, nothing is refitted.

Outputs: results/km_sensitivity/{km_sensitivity.csv, provenance.json}

Usage (from analysis/):  python scripts/km_sensitivity.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from fbaval import data, simulate
from fbaval.provenance import write_provenance

from fbaval.paths import CURATED, DRAFT, RESULTS

OUT = RESULTS / "km_sensitivity"
MODELS = {"draft": DRAFT, "curated": CURATED}
KMS = (0.025, 0.26, 4.5)
DT, T_END = 0.5, 150.0
THR = simulate.INPUTS["diauxie_threshold_mM"][0]


def t50(sim: pd.DataFrame) -> float:
    return float(np.interp(0.5 * sim.biomass_gl.max(), sim.biomass_gl.cummax(), sim.time_h))


def main() -> None:
    fits = pd.read_csv(RESULTS / "flask_dynamics" / "metrics.csv")
    fits = fits[fits.setting == "lag_each"]
    rows = []
    for key, fname in MODELS.items():
        m = simulate.load(str(fname), normalize=True)
        vmax = simulate.vmax_from_mu_max(m, simulate.INPUTS["mu_max_per_h"][0])
        for cond, spec in data.VINESTOCK.items():
            sugars = {s: spec[f"{s}_gl"] for s in ("glucose", "fructose") if spec[f"{s}_gl"] > 0}
            x0 = float(fits[(fits.model == key) & (fits.condition == cond)].x0_gl.iloc[0])
            for km in KMS:
                sim = simulate.run_dfba(m, sugars, x0, T_END, DT, vmax, km, data.VINESTOCK_OTR_MAX_MMOL_L_H, diauxie_threshold_mM=THR)
                rows.append({"model": key, "condition": cond, "km_mM": km, "x0_gl": x0, "x_final_gl": float(sim.biomass_gl.max()), "t50_h": t50(sim)})
    df = pd.DataFrame(rows)
    ref = df[df.km_mM == 0.26].set_index(["model", "condition"])
    df["dx_final_pct"] = [100.0 * (r.x_final_gl / ref.loc[(r.model, r.condition), "x_final_gl"] - 1.0) for r in df.itertuples()]
    df["dt50_h"] = [r.t50_h - ref.loc[(r.model, r.condition), "t50_h"] for r in df.itertuples()]
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "km_sensitivity.csv", index=False)
    write_provenance(OUT, params={"engine": "fbaval.simulate.run_dfba (gem/dfba.py)", "models": {k: v.name for k, v in MODELS.items()}, "km_mM": list(KMS), "dt_h": DT, "t_end_h": T_END,
                                  "x0": "per-condition fit, results/flask_dynamics/metrics.csv", "diauxie_threshold_mM": THR, "fitted_parameters": "none"})
    print(df.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
