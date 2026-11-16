"""Why does the Curated model overpredict final biomass in the batch fermentations? Two tests on existing data and the model.

1. Oxygen balance. Measured oxygen uptake (online OUR = kLa (C* - C)) is integrated over each fermentation and compared
   with what three explanations imply for the missing biomass:
     respired      the carbon left as CO2, so more O2 than predicted was used
     products      the carbon stayed unoxidised in unmeasured products, so less O2 than predicted was used
     not recovered the biomass formed as predicted but was lost before weighing, so O2 is as predicted
   Electron balance: sugar electrons = biomass electrons + 4 O2. Electrons per g biomass come from the model's own end state.
2. Energy yield. The P/O ratio of the model's respiratory chain is lowered by raising H+ per ATP at ATP synthase, and the aerobic
   biomass yield is recomputed, to find the P/O at which the model reaches the measured yields.

Nothing is fitted to the time courses. Outputs: results/overprediction/{oxygen_balance.csv, do_balance.csv, po_sensitivity.csv, po_dfba.csv, provenance.json}

Usage (from analysis/):  python scripts/overprediction_tests.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from fbaval import data, simulate
from fbaval.paths import CURATED, RESULTS
from fbaval.provenance import write_provenance

OUT = RESULTS / "overprediction"
MODEL = CURATED
DT, T_END = 0.5, 150.0
MW_O2 = 32.0
ELECTRONS_PER_HEXOSE = 24.0
BASELINE_H = 15.0                # no measurable growth before this time, so OUR should be about zero
TARGET_PO = (None, 2.25, 2.0, 1.75, 1.5)
CARRIER_H = 1.0                  # r251 moves one proton per ATP exported from the mitochondrion


def integrate_our(cond: str) -> pd.DataFrame:
    """Total O2 used per replicate, as recorded and after removing each sensor's pre-growth offset."""
    rows = []
    for flask, d in data.load_vinestock_presens(cond).groupby("flask"):
        d = d.dropna(subset=["our_g_l_h"]).sort_values("time_h")
        our = d.our_g_l_h.to_numpy() / MW_O2 * 1000.0          # mmol O2 L-1 h-1
        base = float(np.median(our[d.time_h.to_numpy() < BASELINE_H]))
        rows.append({"condition": cond, "flask": flask, "our_peak_mmol_l_h": float(our.max()), "our_baseline_mmol_l_h": base,
                     "o2_raw_mmol_l": float(np.trapezoid(our, d.time_h)), "o2_corrected_mmol_l": float(np.trapezoid(our - base, d.time_h))})
    return pd.DataFrame(rows)


def oxygen_balance(m, x0s: dict[str, float]) -> pd.DataFrame:
    vmax = simulate.vmax_from_mu_max(m, simulate.INPUTS["mu_max_per_h"][0])
    thr = simulate.INPUTS["diauxie_threshold_mM"][0]
    rows = []
    for cond in data.PRIMARY_VINESTOCK:
        spec = data.VINESTOCK[cond]
        sugars = {s: spec[f"{s}_gl"] for s in ("glucose", "fructose") if spec[f"{s}_gl"] > 0}
        sim = simulate.run_dfba(m, sugars, x0s[cond], T_END, DT, vmax, simulate.INPUTS["km_sugar_mM"][0],
                                data.VINESTOCK_OTR_MAX_MMOL_L_H, diauxie_threshold_mM=thr)
        x_pred = float(sim.biomass_gl.max())
        o2_pred = abs(float(sim.o2_mmol_L.iloc[-1]))
        hexose = sum(c / simulate.MW[s] * 1000.0 for s, c in sugars.items())
        electrons = ELECTRONS_PER_HEXOSE * hexose
        e_per_g = (electrons - 4.0 * o2_pred) / x_pred
        fin = data.load_vinestock_final(cond).biomass_gl
        x_meas = float(fin.mean())
        our = integrate_our(cond)
        rows.append({"condition": cond, "x_pred_gl": x_pred, "x_meas_gl": x_meas, "o2_model_mmol_l": o2_pred, "electrons_per_g_biomass": e_per_g,
                     "o2_if_respired": (electrons - x_meas * e_per_g) / 4.0, "o2_if_products": o2_pred * x_meas / x_pred, "o2_if_not_recovered": o2_pred,
                     "o2_measured_raw_mean": our.o2_raw_mmol_l.mean(), "o2_measured_raw_sd": our.o2_raw_mmol_l.std(),
                     "o2_measured_corrected_mean": our.o2_corrected_mmol_l.mean(), "o2_measured_corrected_sd": our.o2_corrected_mmol_l.std(),
                     "our_peak_measured": our.our_peak_mmol_l_h.mean(), "our_baseline": our.our_baseline_mmol_l_h.mean()})
    return pd.DataFrame(rows)


# F1F0-ATPase and its reverse reaction, which share one H+/ATP. r251 is the ADP/ATP translocase, not the synthase.
SYNTHASE = ("r245", "r246")


def yield_at(m, sugar: str, q: float) -> float:
    with m:
        simulate.set_sugars(m, {sugar: q})
        mu = m.slim_optimize(error_value=0.0) * simulate.biomass_weight(m)
    return mu / (q * simulate.MW[sugar] / 1000.0)


def max_atp(m) -> float:
    import curate_model as cm
    with m:
        cm.set_medium(m, cm.MINIMAL_GLUCOSE)
        return float(cm.max_atp_yield(m))


def po_sensitivity(m) -> pd.DataFrame:
    """P/O for NADH = protons pumped per NADH (complexes I, III, IV) / protons per cytosolic ATP (synthase + carrier)."""
    syn = m.reactions.get_by_id(SYNTHASE[0])
    h_po = m.metabolites.H_PO
    h_base = abs(syn.metabolites[h_po])
    pumped = sum(abs(m.reactions.get_by_id(r).metabolites[h_po]) for r in ("r239", "r243", "r244"))
    vmax = simulate.vmax_from_mu_max(m, simulate.INPUTS["mu_max_per_h"][0])
    po_base = pumped / (h_base + CARRIER_H)
    rows = []
    for target in TARGET_PO:
        po = po_base if target is None else target
        mod = m if target is None else scaled_model(m, po)
        vm = simulate.vmax_from_mu_max(mod, simulate.INPUTS["mu_max_per_h"][0])
        rows.append({"po_nadh": po, "po_succinate": (pumped - 4.0) * po / pumped, "h_per_cytosolic_atp": pumped / po, "max_atp_per_glucose": max_atp(mod),
                     "yield_glucose_gg": yield_at(mod, "glucose", vm), "yield_fructose_gg": yield_at(mod, "fructose", vm)})
    return pd.DataFrame(rows)


def scaled_model(m, po: float):
    """A copy of the model with the F1F0-ATPase pair set to the given P/O for NADH."""
    h_po = m.metabolites.H_PO
    h_base = abs(m.reactions.get_by_id(SYNTHASE[0]).metabolites[h_po])
    pumped = sum(abs(m.reactions.get_by_id(r).metabolites[h_po]) for r in ("r239", "r243", "r244"))
    h_new = pumped / po - CARRIER_H
    out = m.copy()
    for rid in SYNTHASE:
        r = out.reactions.get_by_id(rid)
        r.add_metabolites({out.metabolites.H_PO: np.sign(r.metabolites[out.metabolites.H_PO]) * (h_new - h_base),
                           out.metabolites.H_POm: np.sign(r.metabolites[out.metabolites.H_POm]) * (h_new - h_base)})
    return out


def po_dfba(m, x0s: dict[str, float]) -> pd.DataFrame:
    """Final biomass in the batch fermentations at a lower P/O, sequential uptake, nothing else changed."""
    thr = simulate.INPUTS["diauxie_threshold_mM"][0]
    rows = []
    for po in (None, 2.25, 2.0):
        mod = m if po is None else scaled_model(m, po)
        vmax = simulate.vmax_from_mu_max(mod, simulate.INPUTS["mu_max_per_h"][0])
        for cond in data.PRIMARY_VINESTOCK:
            spec = data.VINESTOCK[cond]
            sugars = {s: spec[f"{s}_gl"] for s in ("glucose", "fructose") if spec[f"{s}_gl"] > 0}
            sim = simulate.run_dfba(mod, sugars, x0s[cond], T_END, DT, vmax, simulate.INPUTS["km_sugar_mM"][0],
                                    data.VINESTOCK_OTR_MAX_MMOL_L_H, diauxie_threshold_mM=thr)
            meas = float(data.load_vinestock_final(cond).biomass_gl.mean())
            x = float(sim.biomass_gl.max())
            rows.append({"po_nadh": po or 12.0 / 4.88, "condition": cond, "x_final_gl": x, "x_meas_gl": meas, "err_pct": 100 * (x / meas - 1)})
    return pd.DataFrame(rows)


def do_balance() -> pd.DataFrame:
    """Total O2 from dissolved O2 and kLa directly: OUR = kLa C*/200 (DO0 - DO), DO in % air saturation, 200 % at equilibrium with the headspace."""
    per_pct = data.VINESTOCK_OTR_MAX_MMOL_L_H / 200.0          # kLa C* / 200, mmol L-1 h-1 per % air saturation
    rows = []
    for cond in data.PRIMARY_VINESTOCK:
        for flask, d in data.load_vinestock_presens(cond).groupby("flask"):
            d = d.dropna(subset=["do_pct_as"]).sort_values("time_h")
            t, do = d.time_h.to_numpy(), d.do_pct_as.to_numpy()
            do0 = float(np.median(do[t < BASELINE_H]))
            rows.append({"condition": cond, "flask": flask, "do0_pct": do0, "do_min_pct": float(do.min()), "o2_from_do_mmol_l": float(np.trapezoid(per_pct * (do0 - do), t)),
                         "o2_from_do_nominal200_mmol_l": float(np.trapezoid(per_pct * (200.0 - do), t))})
    return pd.DataFrame(rows)


def main() -> None:
    m = simulate.load(str(MODEL), normalize=True)
    fits = pd.read_csv(RESULTS / "flask_dynamics" / "metrics.csv")
    fits = fits[(fits.setting == "lag_each") & (fits.model == "curated")].set_index("condition").x0_gl.to_dict()
    bal = oxygen_balance(m, fits)
    po = po_sensitivity(m)
    dfba = po_dfba(m, fits)
    dob = do_balance()
    OUT.mkdir(parents=True, exist_ok=True)
    dfba.to_csv(OUT / "po_dfba.csv", index=False)
    dob.to_csv(OUT / "do_balance.csv", index=False)
    bal.to_csv(OUT / "oxygen_balance.csv", index=False)
    po.to_csv(OUT / "po_sensitivity.csv", index=False)
    write_provenance(OUT, params={"engine": "fbaval.simulate.run_dfba (gem/dfba.py)", "model": MODEL.name, "dt_h": DT, "x0": "per-condition fit, results/flask_dynamics/metrics.csv",
                                  "our": "online our_g_l_h, kLa 96 1/h", "baseline_h": BASELINE_H, "target_po": [p for p in TARGET_PO if p], "fitted_parameters": "none"})
    pd.set_option("display.width", 220, "display.max_columns", 30)
    print(bal.round(2).to_string(index=False))
    print(po.round(3).to_string(index=False))
    print(dfba.round(2).to_string(index=False))
    print(dfba.assign(a=dfba.err_pct.abs()).groupby("po_nadh").a.mean().round(1).to_string())
    print(dob.round(1).to_string(index=False))
    print(dob.groupby("condition", sort=False)[["o2_from_do_mmol_l"]].agg(["mean", "std"]).round(0).to_string())


if __name__ == "__main__":
    main()
