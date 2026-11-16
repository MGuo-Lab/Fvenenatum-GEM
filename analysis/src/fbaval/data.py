"""Loaders for the batch fermentation data (Vinestock & Guo 2026).

Shake flasks of F. venenatum A3/5: 50 mL RHM medium in 250 mL baffled flasks, 28 °C, 220 rpm, humidified headspace
with 200 % of the O2 of air, kLa 96 1/h, three flasks (A1-A3) per condition. Each condition has a folder in
``data/flask`` with the HPLC and end-point workbook (``workbook.xlsx``) and the processed PreSens online series
(``Processed_Presens_Data.xlsx``, with ``Processed_HPLC_Data.xlsx`` for the time of the first HPLC sample).
"""

from __future__ import annotations

import re
from pathlib import Path

import openpyxl
import pandas as pd

from fbaval.paths import DATA_DIR, FLASK_DIR  # noqa: F401  (DATA_DIR is part of the public interface)

# Maximum oxygen transfer rate of the shake flasks: 1.50 gO2/L/h at kLa 96 1/h with a 200 %-air headspace
# -> 46.9 mmol/L/h.
VINESTOCK_KLA_PER_H = 96.0
VINESTOCK_OTR_MAX_MMOL_L_H = 1.50 / 32.0 * 1000.0

# presens_dir: name of the run folder of the raw PreSens exports (FBAVAL_RAW_DIR/presens/<presens_dir>)
VINESTOCK = {
    "Glc15": {"presens_dir": "19012026Glucose15gpL", "glucose_gl": 15.0, "fructose_gl": 0.0},
    "Glc10Fru5": {"presens_dir": "02122025GluFru10gpl5gpl", "glucose_gl": 10.0, "fructose_gl": 5.0},
    "Glc5Fru10": {"presens_dir": "15122025FruGlu10gpl5gpl", "glucose_gl": 5.0, "fructose_gl": 10.0},
    "Fru15": {"presens_dir": "12012026Fructose15gpL", "glucose_gl": 0.0, "fructose_gl": 15.0},
}

# The four conditions, in the order used in every table and figure.
PRIMARY_VINESTOCK = ("Glc15", "Glc10Fru5", "Glc5Fru10", "Fru15")

_CONC = re.compile(r"^(?P<compound>.+?) Conc \(g/l\)")


def condition_dir(condition: str, data_dir: Path | None = None) -> Path:
    """Folder of one condition, ``data/flask/<condition>``."""
    return Path(data_dir or FLASK_DIR) / condition


def _workbook(condition: str, data_dir: Path | None):
    return openpyxl.load_workbook(condition_dir(condition, data_dir) / "workbook.xlsx", data_only=True)


def load_vinestock_hplc(condition: str, data_dir: Path | None = None) -> pd.DataFrame:
    """Long table: time_h, flask (1-3), compound (lower case), conc_gl.

    Flasks are assigned by column order within each compound, because some
    header labels repeat 'Flask A1' for all three flasks.
    """
    rows = list(_workbook(condition, data_dir)["HPLC"].iter_rows(values_only=True))
    header = rows[0]
    t_col = header.index("Elapsed Time (h)")
    groups: dict[str, list[int]] = {}
    for i, h in enumerate(header):
        m = _CONC.match(str(h or ""))
        if m:
            groups.setdefault(m.group("compound").strip().lower(), []).append(i)
    out = []
    for row in rows[1:]:
        if not (isinstance(row[0], str) and row[0].startswith("TP")):
            continue
        for compound, cols in groups.items():
            for flask, i in enumerate(cols, start=1):
                v = row[i]
                if isinstance(v, (int, float)):
                    out.append((float(row[t_col]), flask, compound, float(v)))
    return pd.DataFrame(out, columns=["time_h", "flask", "compound", "conc_gl"])


_PRODUCT_DIVISOR = {"glycerol": "Glycerol Calibration Divisor", "ethanol": "Ethanol Calibration Divisor",
                    "acetate": "Acetate Calibration Divisor", "lactate": "Lactic acid Calibration Divisor"}


def load_vinestock_products(condition: str, data_dir: Path | None = None) -> pd.DataFrame:
    """HPLC products recomputed from peak areas and each compound's own calibration divisor.

    Long table: time_h, flask (1-3), compound, area, conc_gl, conc_stored_gl (the workbook value).
    ``conc_gl`` always uses the divisor named for the compound.

    Do not read "lactate" as lactic acid. The workbook lactate areas are a peak at 16.1 min, and lactic
    acid elutes mainly at 19.1 min, where the fermentations show no signal (analysis/scripts/flask_hplc.py).
    """
    wb = _workbook(condition, data_dir)
    divisor = {r[0]: r[1] for r in wb["Calibration"].iter_rows(values_only=True) if r and r[0]}
    rows = list(wb["HPLC"].iter_rows(values_only=True))
    header = [str(h or "") for h in rows[0]]
    t_col = header.index("Elapsed Time (h)")
    num = lambda v: float(v) if isinstance(v, (int, float)) else 0.0
    out = []
    for compound, label in _PRODUCT_DIVISOR.items():
        name = compound.capitalize()
        areas = [j for j, h in enumerate(header) if h.startswith(f"{name} Peak Area")]
        concs = [j for j, h in enumerate(header) if h.startswith(f"{name} Conc")]
        for row in rows[1:]:
            if not (isinstance(row[0], str) and row[0].startswith("TP")):
                continue
            for flask, (a, c) in enumerate(zip(areas, concs), start=1):
                out.append((float(row[t_col]), flask, compound, num(row[a]), num(row[a]) / float(divisor[label]), num(row[c])))
    return pd.DataFrame(out, columns=["time_h", "flask", "compound", "area", "conc_gl", "conc_stored_gl"])


def load_vinestock_final(condition: str, data_dir: Path | None = None) -> pd.DataFrame:
    """Per-flask end-point values: biomass_gl, yield_gg, end_volume_ml, end_ph."""
    ws = _workbook(condition, data_dir)["Final Data"]
    table = {row[0]: row[1:4] for row in ws.iter_rows(values_only=True) if row and row[0]}
    keys = {"biomass_gl": "Biomass Concentration (mg/ml = g/l)",
            "yield_gg": "Biomass Yield (g/g)",
            "end_volume_ml": "End Broth Volume (ml)",
            "end_ph": "End pH"}
    df = pd.DataFrame({k: [float(str(v).rstrip("%")) / (100 if "%" in str(v) else 1)
                           for v in table[label]] for k, label in keys.items()})
    df.insert(0, "flask", [1, 2, 3])
    return df


_PRESENS_COLS = {
    "biomass_gl": re.compile(r"^Biomass Concentration \(g/l\) \(Position (\d)\)$"),
    "ph": re.compile(r"^pH \(Position (\d)\)$"),
    "do_pct_as": re.compile(r"^O2 .*Position (\d)"),
    "our_g_l_h": re.compile(r"^OUR \(g/l/h\) \(Position (\d)\)$"),
}


def load_vinestock_presens(condition: str, data_dir: Path | None = None) -> pd.DataFrame:
    """PreSens SFR vario online series, one row per reading and flask: time_h,
    flask (1-3), biomass_gl, ph, do_pct_as (dissolved O2, % air saturation),
    our_g_l_h. Biomass and pH are calibrated to each flask's start and end-point values.

    time_h counts from the first HPLC sample, as in the HPLC tables. It comes from the absolute timestamps,
    because the file's "Time Lapsed (h)" starts at the first online reading, up to 1 h after that sample."""
    folder = condition_dir(condition, data_dir)
    ws = openpyxl.load_workbook(folder / "Processed_Presens_Data.xlsx", read_only=True,
                                data_only=True).worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    header = [str(h or "") for h in rows[0]]
    t_col = header.index("Timestamp, absolute")
    hplc = list(openpyxl.load_workbook(folder / "Processed_HPLC_Data.xlsx", read_only=True,
                                       data_only=True).worksheets[0].iter_rows(values_only=True))
    t0 = hplc[1][list(hplc[0]).index("Date Time")]
    cols: dict[tuple[int, str], int] = {}
    for i, h in enumerate(header):
        for var, pat in _PRESENS_COLS.items():
            m = pat.match(h)
            if m:
                cols[(int(m.group(1)), var)] = i
    num = lambda v: float(v) if isinstance(v, (int, float)) else float("nan")
    out = []
    for row in rows[1:]:
        row = tuple(row) + (None,) * (len(header) - len(row))   # read-only rows drop trailing blanks
        if row[t_col] is None:
            continue
        for flask in (1, 2, 3):
            out.append([(row[t_col] - t0).total_seconds() / 3600.0, flask] +
                       [num(row[cols[(flask, var)]]) for var in _PRESENS_COLS])
    return pd.DataFrame(out, columns=["time_h", "flask", *_PRESENS_COLS])
