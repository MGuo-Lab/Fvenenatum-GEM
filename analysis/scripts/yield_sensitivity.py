"""Sensitivity of the aerobic biomass yield to biomass composition and maintenance, Draft and Curated.

Each macromolecular mass fraction (BIOMASS_COMPOSITION in gem/curate_model.py) is changed by +/-25 % with the others rescaled so
that biomass stays 1 g/mmol. Growth- and non-growth-associated ATP costs are changed by +/-25 %. Yield is
computed at the sugar uptake that gives mu_max 0.20 1/h on the unchanged model, with O2 not limiting.
For each ATP cost the script also solves for the change that would reach the measured yield.

Outputs: results/yield_sensitivity/{sensitivity.csv, closure.csv, provenance.json}

Usage (from analysis/):  python scripts/yield_sensitivity.py
"""

from __future__ import annotations

import curate_model as cm
import pandas as pd
from scipy.optimize import brentq

from fbaval import data, simulate
from fbaval.provenance import write_provenance

from fbaval.paths import CURATED, DRAFT, RESULTS

OUT = RESULTS / "yield_sensitivity"
MODELS = {"draft": DRAFT, "curated": CURATED}
STEP = 0.25
FIX = {"Protein": "Protein", "CELLWALL": "Cell wall", "LIPIDS": "Lipids", "RNA": "RNA", "DNA": "DNA", "POOL": "Small-molecule pool"}


def yield_gg(m, sugar: str, q: float) -> float:
    with m:
        simulate.set_sugars(m, {sugar: q})
        mu = m.slim_optimize(error_value=0.0) * simulate.biomass_weight(m)
    return mu / (q * simulate.MW[sugar] / 1000.0)


def scale_fraction(m, comp: str, factor: float) -> None:
    """Scale one mass fraction, rescaling the others so the fractions still sum to one."""
    total = sum(cm.BIOMASS_COMPOSITION.values())
    f = {k: v / total for k, v in cm.BIOMASS_COMPOSITION.items()}
    new_i = f[comp] * factor
    rest = (1.0 - new_i) / (1.0 - f[comp])
    g = m.reactions.GROWTH
    for k in f:
        met = m.metabolites.get_by_id(k)
        ratio = factor if k == comp else rest
        g.add_metabolites({met: g.metabolites[met] * (ratio - 1.0)})


def scale_gam(m, factor: float) -> None:
    g = m.reactions.GROWTH
    for mid in ("ATP", "H2O", "ADP", "PI"):
        met = m.metabolites.get_by_id(mid)
        g.add_metabolites({met: g.metabolites[met] * (factor - 1.0)})


def scale_ngam(m, factor: float) -> None:
    r = m.reactions.NGAM
    base = r.lower_bound
    r.upper_bound = max(r.upper_bound, base * factor)
    r.lower_bound = base * factor
    r.upper_bound = base * factor


def main() -> None:
    measured = {s: float(data.load_vinestock_final(c).biomass_gl.mean()) / 15.0 for s, c in (("glucose", "Glc15"), ("fructose", "Fru15"))}
    rows, closure = [], []
    for key, fname in MODELS.items():
        base = simulate.load(str(fname), normalize=True)
        q = simulate.vmax_from_mu_max(base, simulate.INPUTS["mu_max_per_h"][0])
        y0 = {s: yield_gg(base, s, q) for s in measured}
        cases = [(f"{v} fraction", lambda m, fac, k=k: scale_fraction(m, k, fac)) for k, v in FIX.items()]
        cases += [("Growth-associated ATP", scale_gam), ("Non-growth-associated ATP", scale_ngam)]
        for label, fn in cases:
            for sign, fac in (("-25%", 1 - STEP), ("+25%", 1 + STEP)):
                m = base.copy()
                fn(m, fac)
                for s in measured:
                    y = yield_gg(m, s, q)
                    rows.append({"model": key, "parameter": label, "change": sign, "sugar": s, "uptake_mmol_g_h": q, "yield_gg": y,
                                 "yield_base_gg": y0[s], "yield_change_pct": 100 * (y / y0[s] - 1),
                                 "yield_measured_gg": measured[s], "error_vs_measured_pct": 100 * (y / measured[s] - 1)})
        for label, fn in (("Growth-associated ATP", scale_gam), ("Non-growth-associated ATP", scale_ngam)):
            for s, target in measured.items():
                def gap(fac):
                    m = base.copy()
                    fn(m, fac)
                    return yield_gg(m, s, q) - target
                try:
                    fac = brentq(gap, 1.0, 6.0, xtol=1e-3)
                except ValueError:
                    fac = float("nan")
                closure.append({"model": key, "parameter": label, "sugar": s, "yield_base_gg": y0[s], "yield_measured_gg": target,
                                "factor_to_reach_measured": fac})
    OUT.mkdir(parents=True, exist_ok=True)
    df, cl = pd.DataFrame(rows), pd.DataFrame(closure)
    df.to_csv(OUT / "sensitivity.csv", index=False)
    cl.to_csv(OUT / "closure.csv", index=False)
    write_provenance(OUT, params={"engine": "cobra static FBA", "models": {k: v.name for k, v in MODELS.items()}, "step": STEP, "fractions": cm.BIOMASS_COMPOSITION,
                                  "mu_max": simulate.INPUTS["mu_max_per_h"][0], "o2": "not limiting", "fitted_parameters": "none"})
    print(df[df.sugar == "glucose"].pivot_table(index=["model", "parameter"], columns="change", values="yield_change_pct").round(2).to_string())
    print(cl.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
