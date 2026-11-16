"""How much of the Curated model's gain in final biomass comes from the rules that set the energy yield.

The Draft model with only these rules applied (reaction directionality and the GTP/GDP carrier direction), run like the Curated
and Draft models in the batch fermentations: sequential uptake, assumed initial biomass, nothing fitted.

Outputs: results/energy_rules_effect/{effects.csv, provenance.json}

Usage (from analysis/):  python scripts/energy_rules_effect.py
"""

from __future__ import annotations

import curate_model as cm
import pandas as pd
from curation_stage_effects import transfer

from fbaval import data, simulate
from fbaval.provenance import write_provenance

from fbaval.paths import CURATED, CURATION_DIR, DRAFT, RESULTS

OUT = RESULTS / "energy_rules_effect"
X0, DT, T_END = 0.01, 0.5, 150.0
VARIANTS = {"Draft": [], "Draft + stage 4 directions + carrier": ["4_thermo", "6_v3_gtp_carrier"],
            "Draft + stage 6 directions + carrier": ["6_v2_thermo", "6_v3_gtp_carrier"],
            "Draft + both direction sets + carrier": ["4_thermo", "6_v2_thermo", "6_v3_gtp_carrier"], "Curated": None}


def run(mod, label: str) -> list[dict]:
    with mod:
        cm.set_medium(mod, cm.MINIMAL_GLUCOSE)
        atp = float(cm.max_atp_yield(mod))
    vmax = simulate.vmax_from_mu_max(mod, simulate.INPUTS["mu_max_per_h"][0])
    rows = []
    for cond in data.PRIMARY_VINESTOCK:
        spec = data.VINESTOCK[cond]
        sugars = {s: spec[f"{s}_gl"] for s in ("glucose", "fructose") if spec[f"{s}_gl"] > 0}
        sim = simulate.run_dfba(mod, sugars, X0, T_END, DT, vmax, simulate.INPUTS["km_sugar_mM"][0], data.VINESTOCK_OTR_MAX_MMOL_L_H,
                                diauxie_threshold_mM=simulate.INPUTS["diauxie_threshold_mM"][0])
        meas = float(data.load_vinestock_final(cond).biomass_gl.mean())
        x = float(sim.biomass_gl.max())
        rows.append({"variant": label, "condition": cond, "max_atp_per_glucose": atp, "x_final_gl": x, "x_meas_gl": meas, "err_pct": 100 * (x / meas - 1)})
    return rows


def main() -> None:
    draft = simulate.load(str(DRAFT), normalize=True)
    cur = simulate.load(str(CURATED), normalize=True)
    log = pd.read_csv(CURATION_DIR / "curation_log.csv")
    rows = []
    for label, stages in VARIANTS.items():
        if stages is None:
            rows += run(cur, label)
            continue
        m = draft.copy()
        transfer(m, cur, sorted(set(log[log.stage.isin(stages)].reaction_id)))
        rows += run(m, label)
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "effects.csv", index=False)
    write_provenance(OUT, params={"engine": "fbaval.simulate.run_dfba (gem/dfba.py)", "draft": DRAFT.name, "curated": CURATED.name, "x0_gl": X0, "dt_h": DT, "uptake": "sequential", "fitted_parameters": "none"})
    print(df.round(2).to_string(index=False))
    print(df.assign(a=df.err_pct.abs()).groupby("variant", sort=False).agg(atp=("max_atp_per_glucose", "first"), mean_abs_err=("a", "mean")).round(2).to_string())


if __name__ == "__main__":
    main()
