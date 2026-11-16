"""Tidy CSV copies of the batch fermentation data, written next to the workbooks in data/flask/.

  hplc.csv     condition, flask, time_h, compound, conc_gl: sugars and products by HPLC. Glycerol, ethanol and acetate
               are recomputed from peak areas with each compound's own calibration (fbaval.data.load_vinestock_products).
               The workbook column "lactate" is left out: it holds a peak at 16.1 min, tentatively assigned to arabitol
               and not lactic acid, whose re-integrated area is in results/flask_hplc/unidentified.csv.
  online.csv   condition, flask, time_h, biomass_gl, ph, do_pct_as, our_g_l_h: PreSens online series. Biomass and pH
               are calibrated to the start and end points of each flask; DO is in % of air saturation.
  final.csv    condition, flask, biomass_gl, yield_gg, end_volume_ml, end_ph: end-point dry weight and yield.

Time is in hours from the first HPLC sample of each fermentation, for HPLC and online data alike.

Usage (from analysis/):  python scripts/export_flask_data.py
"""

from __future__ import annotations

import pandas as pd

from fbaval import data
from fbaval.paths import FLASK_DIR

SUGARS = ("glucose", "fructose")
PRODUCTS = ("glycerol", "ethanol", "acetate")


def main() -> None:
    hplc, online, final = [], [], []
    for cond in data.PRIMARY_VINESTOCK:
        sugars = data.load_vinestock_hplc(cond)
        sugars = sugars[sugars.compound.isin(SUGARS)][["time_h", "flask", "compound", "conc_gl"]]
        prod = data.load_vinestock_products(cond)
        prod = prod[prod.compound.isin(PRODUCTS)][["time_h", "flask", "compound", "conc_gl"]]
        hplc.append(pd.concat([sugars, prod]).assign(condition=cond))
        online.append(data.load_vinestock_presens(cond).assign(condition=cond))
        final.append(data.load_vinestock_final(cond).assign(condition=cond))
    hplc = pd.concat(hplc).sort_values(["condition", "flask", "compound", "time_h"], key=lambda s: s.map(
        {c: i for i, c in enumerate(data.PRIMARY_VINESTOCK)}) if s.name == "condition" else s)
    online = pd.concat(online)
    final = pd.concat(final)
    hplc[["condition", "flask", "time_h", "compound", "conc_gl"]].round(4).to_csv(FLASK_DIR / "hplc.csv", index=False)
    online[["condition", "flask", "time_h", "biomass_gl", "ph", "do_pct_as", "our_g_l_h"]].round(4).to_csv(FLASK_DIR / "online.csv", index=False)
    final[["condition", "flask", "biomass_gl", "yield_gg", "end_volume_ml", "end_ph"]].round(4).to_csv(FLASK_DIR / "final.csv", index=False)
    print(f"hplc {len(hplc)} rows, online {len(online)} rows, final {len(final)} rows -> {FLASK_DIR}")


if __name__ == "__main__":
    main()
