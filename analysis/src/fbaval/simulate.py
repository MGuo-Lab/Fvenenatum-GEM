"""Model-vs-data simulations without parameter fitting.

All inputs are fixed before any comparison and labelled by origin:
``measured`` (from the experiment description), ``literature``, or
``assumed``. Nothing here is tuned to the outputs being compared.

Static FBA gives structural predictions (yield, which byproducts appear
under O2 limitation). dFBA uses ``run_dfba_multi`` (gem/dfba.py) with Michaelis-
Menten uptake and a volumetric O2 transfer limit, qO2_max(t) = OTR_max / X(t).
"""

from __future__ import annotations

import cobra
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis, pfba

import curate_model as cm
from dfba import run_dfba_multi

MW = {"glucose": 180.156, "fructose": 180.156, "sucrose": 342.30, "xylose": 150.13,
      "ethanol": 46.07, "glycerol": 92.09, "acetate": 60.05, "pyruvate": 88.06,
      "succinate": 118.09, "co2": 44.01}
EXCHANGE = {"glucose": "BOUNDARY_GLCe", "fructose": "BOUNDARY_FRUe",
            "sucrose": "BOUNDARY_SUCe", "xylose": "BOUNDARY_XYLe",
            "ethanol": "BOUNDARY_ETHe", "glycerol": "BOUNDARY_GLe",
            "acetate": "BOUNDARY_ACe", "pyruvate": "BOUNDARY_PYRe",
            "succinate": "BOUNDARY_SUCCe", "co2": "BOUNDARY_CO2e", "o2": "BOUNDARY_O2e"}
BYPRODUCTS = ["ethanol", "glycerol", "acetate", "pyruvate", "succinate"]
NH4_EXCHANGE = "BOUNDARY_NH4e"
NH4_VMAX, NH4_KM_MM = 20.0, 0.01     # assumed: minimal-medium bound, near-zero Km

# Literature-derived inputs (not fitted to the comparison data).
INPUTS = {
    "mu_max_per_h": (0.20, "literature", "Quorn process 0.17-0.20 1/h, Wiebe 2002"),
    "km_sugar_mM": (0.26, "assumed", "between the glucose (0.025 mM) and fructose (4.5 mM) values of A. niger MstA, van Kuyk 2004"),
    "fructose_uptake_as_glucose": (True, "assumed", "no F. venenatum fructose kinetics published"),
    "vinestock_otr_mmol_L_h": (46.9, "measured", "1.50 gO2/L/h at kLa 96 1/h, 200% air (Vinestock & Guo 2026)"),
    "air_o2_saturation_mM": (0.24, "literature", "O2 solubility in water at 28 C, 21% O2"),
    "diauxie_threshold_mM": (0.5, "assumed", "glucose-first uptake: fructose blocked while glucose > 0.5 mM"),
}


def load(path: str, normalize: bool = False) -> cobra.Model:
    """Read an SBML model; ``normalize`` makes 1 unit of biomass flux = 1 g."""
    model = cobra.io.read_sbml_model(path)
    if normalize:
        normalize_biomass(model)
    return model


def biomass_weight(model) -> float:
    return cm.biomass_weight(model) or 1.0


def normalize_biomass(model, tol: float = 1e-3) -> float:
    """Scale the reaction that makes BIOMASS so one unit of it weighs 1 g.

    A pure unit change: yields in g/g are unchanged, but the biomass flux
    becomes a specific growth rate in 1/h, which dFBA uses directly. The
    curated models are already normalised (weight 1.0) and are left alone.
    Returns the weight before scaling (g/mmol)."""
    bw = biomass_weight(model)
    if abs(bw - 1.0) > tol:
        for r in model.reactions:
            if "BIOMASS" in {m.id for m in r.products}:
                # every input (and GAM) per mmol BIOMASS shrinks to what 1 g needs
                r.add_metabolites({m: c * (1.0 / bw - 1.0) for m, c in r.metabolites.items()
                                   if m.id != "BIOMASS"})
    return bw


def set_sugars(model, uptakes: dict[str, float], o2_cap: float | None = None,
               forced: bool = False) -> None:
    """Minimal medium with the given sugar uptakes (mmol/gDW/h, positive)."""
    medium = dict(cm.MINIMAL_GLUCOSE)
    medium.pop("BOUNDARY_GLCe")
    medium["BOUNDARY_O2e"] = -1000.0 if o2_cap is None else -abs(o2_cap)
    for sugar, q in uptakes.items():
        medium[EXCHANGE[sugar]] = -q
    cm.set_medium(model, medium)
    for sugar, q in uptakes.items():
        rxn = model.reactions.get_by_id(EXCHANGE[sugar])
        rxn.upper_bound = -q if forced else min(rxn.upper_bound, 1000.0)


def static_profile(model, sugar: str, q: float,
                   o2_fractions=(1.0, 0.8, 0.6, 0.4, 0.2, 0.1)) -> pd.DataFrame:
    """Growth, yield and byproducts as the O2 supply is reduced.

    Sugar uptake is fixed at ``q`` (forced), O2 is capped at a fraction of
    the unconstrained demand. Yield uses the model's biomass weight.
    """
    bw = biomass_weight(model)
    with model:
        set_sugars(model, {sugar: q})
        demand = -pfba(model).fluxes["BOUNDARY_O2e"]
    rows = []
    ids = [EXCHANGE[b] for b in BYPRODUCTS if EXCHANGE[b] in model.reactions]
    for frac in o2_fractions:
        with model:
            set_sugars(model, {sugar: q}, o2_cap=frac * demand, forced=True)
            sol = model.optimize()
            row = {"sugar": sugar, "q_sugar": q, "o2_fraction": frac,
                   "qO2": frac * demand}
            if sol.status != "optimal":
                row["status"] = sol.status
                rows.append(row)
                continue
            mu = sol.objective_value * bw
            fva = flux_variability_analysis(model, ids, fraction_of_optimum=1.0)
            row.update({"status": "optimal", "mu_per_h": mu,
                        "yield_gg": mu / (q * MW[sugar] / 1000.0)})
            for b in BYPRODUCTS:
                rid = EXCHANGE[b]
                if rid in ids:
                    row[f"{b}_min"] = fva.loc[rid, "minimum"]
                    row[f"{b}_max"] = fva.loc[rid, "maximum"]
            rows.append(row)
    return pd.DataFrame(rows)


def biomass_carbon_fraction(model, sugar: str = "glucose") -> float:
    """Carbon mass fraction of the model's biomass, from an aerobic pFBA
    carbon balance (all C not leaving through exchanges ends in biomass)."""
    with model:
        set_sugars(model, {sugar: 1.0})
        sol = pfba(model)
        c_out = 0.0
        for r in model.reactions:
            if r.id.startswith("BOUNDARY_") and r.id != EXCHANGE[sugar]:
                c_out += sol.fluxes[r.id] * next(iter(r.metabolites)).elements.get("C", 0)
        c_in = -sol.fluxes[EXCHANGE[sugar]] * next(iter(model.reactions.get_by_id(EXCHANGE[sugar]).metabolites)).elements["C"]
        return (c_in - c_out) * 12.011 / 1000.0 / (sol.fluxes["DRAIN_Biomass"] * biomass_weight(model))


def o2_per_sugar(model, sugar: str = "glucose") -> float:
    """mol O2 consumed per mol sugar in unconstrained aerobic growth (pFBA)."""
    with model:
        set_sugars(model, {sugar: 1.0})
        sol = pfba(model)
        return -sol.fluxes["BOUNDARY_O2e"] / -sol.fluxes[EXCHANGE[sugar]]


def vmax_from_mu_max(model, mu_max: float, sugar: str = "glucose") -> float:
    """Uptake Vmax (mmol/gDW/h) that supports the literature mu_max at the
    model's own aerobic yield; derived from model + literature, not data."""
    q = 1.0
    with model:
        set_sugars(model, {sugar: q})
        mu = model.slim_optimize() * biomass_weight(model)
    return mu_max / mu * q


def run_dfba(model, sugars_gl: dict[str, float], x0_gl: float, t_end: float,
             dt: float, vmax: float, km_mM: float, otr_mmol_L_h: float | None,
             byproducts=BYPRODUCTS, diauxie_threshold_mM: float | None = None,
             qo2_cap: float | None = None, vmax_by_sugar: dict[str, float] | None = None,
             ph_kwargs: dict | None = None, nh4_mM: float | None = None,
             coupling: str = "scale_all") -> pd.DataFrame:
    """dFBA in g/L (volume 1 L). qO2_max = OTR / X when ``otr_mmol_L_h``.

    With ``diauxie_threshold_mM`` and glucose present, every other sugar is
    blocked while glucose is above the threshold (glucose-first uptake).
    ``qo2_cap`` (mmol/gDW/h) instead applies a constant specific O2 cap; it
    ignores biomass growth.
    ``vmax_by_sugar`` overrides ``vmax`` per sugar.
    ``ph_kwargs`` switches on the pH module of gem/dfba.py (keyword arguments of
    ``run_dfba_multi``: pH_params, initial_pH, medium_phosphate_mM, ...);
    secreted acetate, pyruvate and succinate count as acids and NH4/P/S
    uptake is tracked. The output then has a ``pH`` column.
    ``nh4_mM`` makes the medium NH4 a finite pool, taken up on demand
    (column ``nh4_mM``); by default NH4 is unlimited.
    ``coupling`` is the ``pH_coupling`` option of ``run_dfba_multi``: "growth_cap" keeps sugar uptake
    at its kinetic rate when growth is held back, by pH or by running out of
    NH4, and sends the surplus carbon to CO2 or, beyond the O2 limit, to
    byproducts.

    The biomass flux is used as the specific growth rate, so the model's
    BIOMASS should be normalised to 1 g/mmol (``load(..., normalize=True)``).
    """
    with model:
        set_sugars(model, {s: 0.0 for s in sugars_gl})
        subs = []
        for s, c in sugars_gl.items():
            rxn = model.reactions.get_by_id(EXCHANGE[s])
            rxn.upper_bound = 0.0
            subs.append({"name": s, "rxn": rxn, "mw": MW[s],
                         "vmax": (vmax_by_sugar or {}).get(s, vmax),
                         "km": km_mM, "conc_mm": c / MW[s] * 1000.0})
        if nh4_mM is not None:
            # uptake capacity of the minimal medium; non-limiting until the pool is gone
            subs.append({"name": "nh4", "rxn": model.reactions.get_by_id(NH4_EXCHANGE),
                         "mw": 18.04, "vmax": NH4_VMAX, "km": NH4_KM_MM,
                         "conc_mm": nh4_mM, "on_demand": True})
        byp_ids = [EXCHANGE[b] for b in byproducts if EXCHANGE[b] in model.reactions]
        byp_ids += [EXCHANGE["co2"], EXCHANGE["o2"]]
        biomass_rxn = model.reactions.get_by_id("DRAIN_Biomass")
        o2_fn = None if otr_mmol_L_h is None else (
            lambda t, x: otr_mmol_L_h / max(x, 1e-6))
        if qo2_cap is not None:
            o2_fn = lambda t, x: qo2_cap
        primary = (EXCHANGE["glucose"] if diauxie_threshold_mM is not None
                   and "glucose" in sugars_gl and len(sugars_gl) > 1 else None)
        extra = {}
        if ph_kwargs:
            extra = {"acid_keys": {EXCHANGE[a]: a for a in ("acetate", "pyruvate", "succinate")
                                   if EXCHANGE[a] in model.reactions}, **ph_kwargs}
        res = run_dfba_multi(model, biomass_rxn, subs, x0_gl, 1.0, t_end, dt, byp_ids,
                             diauxic_primary_id=primary,
                             diauxic_threshold_mm=diauxie_threshold_mM or 0.5,
                             o2_cap_fn=o2_fn, continue_on_infeasible=True,
                             pH_coupling=coupling, **extra)
    df = pd.DataFrame(res)
    out = pd.DataFrame({"time_h": df["time"], "biomass_gl": df["biomass"]})
    for s in sugars_gl:
        out[f"{s}_gl"] = df[f"c_{EXCHANGE[s]}"] * MW[s] / 1000.0
    for b in byproducts:
        col = f"c_{EXCHANGE[b]}"
        if col in df:
            out[f"{b}_gl"] = df[col] * MW[b] / 1000.0
    for extra in ("co2", "o2"):
        col = f"c_{EXCHANGE[extra]}"
        if col in df:
            out[f"{extra}_mmol_L"] = df[col]
            out[f"q{extra}"] = df[f"v_{EXCHANGE[extra]}"]
    out["infeasible"] = df["infeasible"].fillna(False) if "infeasible" in df else False
    if nh4_mM is not None:
        out["nh4_mM"] = df[f"c_{NH4_EXCHANGE}"]
    if "coupling_fallback" in df:
        out["coupling_fallback"] = df["coupling_fallback"].fillna(False).astype(bool)
    if ph_kwargs:
        out["pH"] = df["pH"]
        for k in ("nh4", "pi", "so4"):
            if f"{k}_assimilated_mM" in df:
                out[f"{k}_assimilated_mM"] = df[f"{k}_assimilated_mM"]
    return out


def interp_rmse(sim: pd.DataFrame, t_obs, y_obs, col: str = "biomass_gl") -> float:
    import numpy as np
    pred = np.interp(t_obs, sim["time_h"], sim[col])
    return float(np.sqrt(np.mean((pred - np.asarray(y_obs)) ** 2)))
