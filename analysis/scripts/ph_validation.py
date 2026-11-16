"""Predict the pH of the batch fermentations without fitting, and compare it
with the online (PreSens) pH series (Fig. 5).

Stage A tests the chemistry alone. It takes the measured biomass (PreSens,
calibrated to dry weight) and the measured organic acid (HPLC acetate; no lactic
acid is detected, see scripts/flask_hplc.py), then predicts pH with versions
v0-v3 of `fbaval.ph`. The N, P and S content of biomass comes from the biomass
reaction of iJL1143.

Outputs: results/ph_validation/{metrics.csv, series.json, provenance.json}

Usage (from analysis/):  python scripts/ph_validation.py
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

import dfba_kinetics as dk
from fbaval import data, ph, simulate
from fbaval.paths import CURATED, RESULTS
from fbaval.provenance import write_provenance

MODEL = CURATED
OUT = RESULTS / "ph_validation"
ACID_MW = {"acetate": 60.05}


def clean(x):
    return [None if not np.isfinite(v) else round(float(v), 4) for v in np.asarray(x, float)]


def flask_inputs(cond: str):
    """Per flask: PreSens time, biomass formed since t0, measured pH, initial pH
    and HPLC acids interpolated onto the PreSens times (mM)."""
    pre = data.load_vinestock_presens(cond)
    hplc = data.load_vinestock_hplc(cond)
    out = {}
    for f in (1, 2, 3):
        d = pre[pre.flask == f].dropna(subset=["ph"]).sort_values("time_h")
        t = d.time_h.to_numpy()
        x = d.biomass_gl.to_numpy()
        formed = np.clip(x - x[0], 0.0, None)
        acids = {}
        for a, mw in ACID_MW.items():
            h = hplc[(hplc.flask == f) & (hplc.compound == a)].sort_values("time_h")
            if len(h):
                acids[a] = np.interp(t, h.time_h, h.conc_gl.clip(lower=0)) / mw * 1000.0
        out[f] = {"t": t, "x": x, "formed": formed, "ph": d.ph.to_numpy(),
                  "ph0": float(d.ph.iloc[0]), "acids": acids}
    return out


def metrics(pred, meas):
    err = pred - meas
    return {"rmse": float(np.sqrt(np.mean(err ** 2))), "end_error": float(err[-1]),
            "max_abs": float(np.max(np.abs(err)))}


def co2_bound(comp_rq: float = 1.16) -> dict:
    """Largest pH effect of respiratory CO2 in the sparged flasks: peak OUR x RQ
    against the 0.2 L/min gas flow per flask sets pCO2; Henry's law and
    pKa1 6.35 give the bicarbonate formed; the phosphate buffer absorbs it."""
    peak_our = max(data.load_vinestock_presens(c).our_g_l_h.max() for c in data.PRIMARY_VINESTOCK)
    co2_mmol_h = peak_our / 32.0 * 1000.0 * comp_rq * 0.050           # mmol/h per 50 mL flask
    gas_mol_h = 0.2 * 60.0 / 24.7                                      # 0.2 L/min at 28 °C
    pco2 = co2_mmol_h / 1000.0 / gas_mol_h                             # atm
    co2_aq = pco2 * 0.029 * 1000.0                                     # mM, kH(28 °C) ~ 0.029 M/atm
    hco3 = co2_aq * 10 ** (6.0 - 6.35)                                 # at the lowest pH reached, 6.0
    beta = 2.303 * ph.PHOSPHATE_MM * (10 ** -6.79 * 1e-6) / (10 ** -6.79 + 1e-6) ** 2
    return {"peak_our_g_l_h": peak_our, "pco2_atm": pco2, "co2_aq_mM": co2_aq,
            "hco3_mM": hco3, "delta_pH": hco3 / beta}


def stage_a_score(comp, version="v3", pka_shift=0.0):
    """Mean RMSE and mean signed end error of one setting over all flasks."""
    old = dk.PKA_TABLE["phosphate"]
    dk.PKA_TABLE["phosphate"] = (old[0], old[1] + pka_shift, old[2])
    try:
        rm, end = [], []
        for cond in data.PRIMARY_VINESTOCK:
            for d in flask_inputs(cond).values():
                p = ph.predict_ph(d["formed"], d["acids"], d["ph0"], ph.VERSIONS[version], comp)
                rm.append(np.sqrt(np.mean((p - d["ph"]) ** 2)))
                end.append(p[-1] - d["ph"][-1])
    finally:
        dk.PKA_TABLE["phosphate"] = old
    return {"rmse": float(np.mean(rm)), "end_bias": float(np.mean(end))}


def checks(model, comp):
    """Effects considered after v3, each computed rather than assumed."""
    out = {"temperature": {
        # phosphate pKa2 temperature coefficient -0.0028 per K (Good et al. 1966), 25 -> 28 °C
        "pka2_shift": -0.0028 * 3, **stage_a_score(comp, pka_shift=-0.0028 * 3)}}
    out["n_content"] = {f"x{s}": stage_a_score(ph.Composition(comp.n * s, comp.p, comp.s))
                        for s in (0.9, 1.0, 1.1)}
    # OUR integrated over time, divided by the model's O2 per g biomass at mu 0.2,
    # should land near the dry weight if OUR could stand in for biomass.
    from cobra.flux_analysis import pfba
    with model:
        simulate.set_sugars(model, {"glucose": 2.0})
        f = pfba(model).fluxes
        o2_per_g = -f["BOUNDARY_O2e"] / model.slim_optimize()
    ratios = {}
    for cond in data.PRIMARY_VINESTOCK:
        pre = data.load_vinestock_presens(cond)
        fin = data.load_vinestock_final(cond).set_index("flask").biomass_gl
        for fl in (1, 2, 3):
            d = pre[pre.flask == fl].sort_values("time_h")
            our = np.nan_to_num(d.our_g_l_h.to_numpy()) / 32.0 * 1000.0
            cum = float(np.sum(0.5 * (our[1:] + our[:-1]) * np.diff(d.time_h.to_numpy())))
            ratios[f"{cond} F{fl}"] = cum / o2_per_g / float(fin[fl])
    out["our_as_biomass"] = {"o2_per_g_model": o2_per_g, "ratio_to_dry_weight": ratios}
    return out


def stage_a(comp):
    rows, series = [], {}
    for cond in data.PRIMARY_VINESTOCK:
        series[cond] = {}
        for f, d in flask_inputs(cond).items():
            preds = {v: ph.predict_ph(d["formed"], d["acids"], d["ph0"], ver, comp)
                     for v, ver in ph.VERSIONS.items()}
            for v, p in preds.items():
                rows.append({"stage": "A", "condition": cond, "flask": f, "version": v,
                             **metrics(p, d["ph"])})
            step = max(1, len(d["t"]) // 60)
            series[cond][f] = {"t": clean(d["t"][::step]), "x": clean(d["x"][::step]),
                               "ph": clean(d["ph"][::step]),
                               **{v: clean(p[::step]) for v, p in preds.items()},
                               "acids_max_mM": {a: float(c.max()) for a, c in d["acids"].items()}}
    return rows, series


def main() -> None:
    model = simulate.load(str(MODEL))
    comp = ph.biomass_composition(model)
    rows_a, series_a = stage_a(comp)
    co2 = co2_bound()
    extra = checks(model, comp)

    OUT.mkdir(parents=True, exist_ok=True)
    met = pd.DataFrame(rows_a)
    met.to_csv(OUT / "metrics.csv", index=False)
    summary = (met[met.stage == "A"].groupby("version")[["rmse", "end_error", "max_abs"]]
               .agg(lambda s: float(np.mean(np.abs(s)))).round(4).to_dict(orient="index"))
    payload = {
        "meta": {"model": MODEL.name,
                 "composition_mmol_per_g": {"N": comp.n, "P": comp.p, "S": comp.s},
                 "medium_mM": {"phosphate": ph.PHOSPHATE_MM, "NH4Cl": ph.NH4CL_MM,
                               "K2SO4": ph.K2SO4_MM, "MgSO4": ph.MGSO4_MM},
                 "ionic_strength_M": ph.ionic_strength(6.6),
                 "pka2_apparent": dk.apparent_pkas(dk.PKA_TABLE["phosphate"],
                                                   ph.ionic_strength(6.6), ph.TEMP_C)[1],
                 "conditions": list(data.PRIMARY_VINESTOCK)},
        "summary_stage_a": summary, "metrics": met.replace({np.nan: None}).to_dict(orient="records"),
        "co2_check": co2, "checks": extra, "stage_a": series_a,
    }
    (OUT / "series.json").write_text(json.dumps(payload, ensure_ascii=False, allow_nan=False))
    write_provenance(OUT, params={"engine": "dfba_kinetics pH chemistry (charge balance)",
                                  "model": MODEL.name,
                                  "versions": {k: vars(v) for k, v in ph.VERSIONS.items()},
                                  "fitted_parameters": "none"})
    print(pd.DataFrame(summary).T)
    print("CO2 bound:", {k: round(v, 4) for k, v in co2.items()})
    print("checks:", json.dumps({k: v for k, v in extra.items() if k != "our_as_biomass"}, indent=1))
    print("OUR/DW:", {k: round(v, 2) for k, v in extra["our_as_biomass"]["ratio_to_dry_weight"].items()})
    print("composition:", comp)


if __name__ == "__main__":
    main()
