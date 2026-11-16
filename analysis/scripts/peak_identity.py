"""Identity of the 16.1 min refractive-index peak of the batch fermentations, from authentic standards.

The peak is compared with standards injected on the same instrument and column under the same protocol
(oven 30 C, 0.001 M H2SO4, 0.4 mL/min), read from the raw LabSolutions files. Four questions:
  1. Retention. Which standard elutes there, once the column drift between May 2025 (standards) and the
     fermentation runs is removed with glucose and glycerol as neutral anchors? The glucose anchor of the May 2025
     sequence is the median glucose retention time of the 1 uL sample injections of that sequence.
  2. UV. RID and photodiode array run in series, so the UV-to-RI ratio of a peak separates carboxylic acids
     from polyols and sugars without any standard at the same concentration.
  3. Competitors. Arabinose has no standard at 30 C. Its position comes from the two temperatures that do exist
     (25 C on this instrument, 80 C in the Phenomenex RHM chart), checked against mannitol and sorbitol.
  4. Size. Approximate concentration, from a single-level response, and only as a range.

Nothing is fitted. The oven temperature and pump pressure of every injection are read from its status log, because
several folders with similar names were run at other temperatures.

Raw files: the Shimadzu LabSolutions data files (.lcd) of the instrument, under FBAVAL_RAW_DIR/lcd (see fbaval.paths);
they are not part of this repository. Reader: fbaval.lcd.
PDA data: 'PDA 3D Raw Data/3D Raw Data' holds one delta-encoded record per scan; values are checked against the
instrument's own Max Plot stream.

Outputs: results/peak_identity/{standards.csv, delay.csv, flask_peaks.csv, anchored_rt.csv, rt_sources.csv,
         arabinose_position.csv, concentration.csv, model_check.csv, fig_traces.csv, fig_uv.csv, summary.json,
         provenance.json}. Supplementary Fig. S3 is drawn from fig_traces.csv and fig_uv.csv by paper_figures.py.

Usage (from analysis/):  python scripts/peak_identity.py
"""

from __future__ import annotations

import json
import re
import struct
from pathlib import Path

import numpy as np
import olefile
import pandas as pd

from fbaval.lcd import _deltas, _tag, _xml, chromatogram, read_injection
from fbaval.paths import CURATED, RAW_DIR, RESULTS
from fbaval.provenance import write_provenance

OUT = RESULTS / "peak_identity"
RAW = RAW_DIR / "lcd"
REP = RAW / "Fermentation_Exp_Repeats_with_CitricAcid"
P1 = REP / "CitricAcid_Calibration_and_Unknown_Peaks_Identification_080525" / "Protocol1_0.001MH2SO4_30C"
LA_SOLO = REP / "Lactic_Acid_Stds_Solo_230525"
MAY_SEQUENCE = REP / "Fru_Xyl_30gl_wCA_10-05-25"      # sample injections of the May 2025 sequence (glucose anchor)
WK16 = RAW / "YH_S_Ethanol" / "YH_TimeSeriesCharacterisation" / "16 WEEK"
FEB24 = RAW / "YH_MB_Sugar_SugarAlcohol_OA"
STD26 = RAW / "TV" / "HPLC Standards" / "20261005_Standards"
FLASK = RAW / "TV" / "Batch Shake Flask Winter 2025-26"
JAN25 = RAW / "Fermentation_Exp_Single_Substrate_Jan25" / "GFX_T0_conc_+xylitol_calibration_10_01_25"

CENTRAL = 0.25   # min each side of the apex: UV/RI is compared on the central part of the peak, the same way for every peak
FLANK = (0.55, 0.85)   # min from the apex: baseline fitted on these two flanks

# Period, compound, folder, file pattern, apex window (min), nominal g/L, acid or not. Delay comes from the acids of that period.
SPEC = [
    ("May 2025", "Arabitol", P1, "Arabitol_1 g*", (15.9, 16.4), 1.0, False),
    ("May 2025", "Arabitol 0.1 g/L", P1, "Arabitol_0.1 g*", (15.9, 16.4), 0.1, False),
    ("May 2025", "Xylitol", P1, "Xylitol_1 g*", (16.4, 16.9), 1.0, False),
    ("May 2025", "Xylitol 0.1 g/L", P1, "Xylitol_0.1 g*", (16.4, 16.9), 0.1, False),
    ("May 2025", "Sorbitol", P1, "Sorbitol_1 g*", (15.0, 15.3), 1.0, False),
    ("May 2025", "Mannitol", P1, "Mannitol_1 g*", (14.6, 15.0), 1.0, False),
    ("May 2025", "Maltitol", P1, "Maltitol_1 g*", (11.5, 12.0), 1.0, False),
    ("May 2025", "Malic acid", P1, "Malic acid_1 g*", (14.2, 14.7), 1.0, True),
    ("May 2025", "Citric acid", P1, "Citric acid_1.4 g*_1_*", (11.7, 12.1), 1.4, True),
    ("May 2025", "Acetic acid", P1, "Acetic acid_1 g*", (23.2, 23.6), 1.0, True),
    ("May 2025", "Lactic acid (mix L1)", P1, "L1_EtOH*", (18.8, 19.2), None, True),
    ("May 2025", "Lactic reagent minor (mix L1)", P1, "L1_EtOH*", (16.1, 16.4), None, True),
    ("May 2025", "Glycerol (mix L1)", P1, "L1_EtOH*", (19.7, 20.1), None, False),
    ("May 2025", "Lactic acid (solo)", LA_SOLO, "LacticAcid_1gl*", (18.7, 19.1), 1.0, True),
    ("Oct 2026", "Glucose", STD26, "20261005_Glu_2gL001*", (13.2, 13.6), 2.0, False),
    ("Oct 2026", "Fructose", STD26, "20261005_Fru_2gL001*", (14.5, 14.9), 2.0, False),
    ("Oct 2026", "Glycerol", STD26, "20261005_Gly_2gL001*", (19.6, 20.0), 2.0, False),
    ("Oct 2026", "Lactic acid", STD26, "20261005_LA_2gL001*", (18.9, 19.3), 2.0, True),
    ("Oct 2026", "Lactic reagent minor 15.4", STD26, "20261005_LA_2gL001*", (15.2, 15.6), None, True),
    ("Oct 2026", "Lactic reagent minor 16.3", STD26, "20261005_LA_2gL001*", (16.1, 16.5), None, True),
    ("Oct 2026", "Acetic acid", STD26, "20261005_Acetic_2gL001*", (23.2, 23.6), 2.0, True),
    ("Dec 2024", "Arabitol", WK16, "*_Arabitol*", (15.8, 16.2), 2.0, False),
    ("Dec 2024", "Arabinose", WK16, "*_Arabinose*", (15.8, 16.2), 10.0, False),
    ("Dec 2024", "Mannitol", WK16, "*_Mannitol*", (14.5, 14.9), 10.0, False),
    ("Dec 2024", "Ribose", WK16, "*_Ribose*", (16.5, 17.0), 10.0, False),
    ("Dec 2024", "Mannose", WK16, "*_Mannose*", (14.1, 14.5), 10.0, False),
    ("Dec 2024", "Galactose", WK16, "*_Galactose*", (14.2, 14.6), 10.0, False),
    ("Dec 2024", "Glucose (mix)", WK16, "*10gL_Glu_Fru*", (13.0, 13.4), 10.0, False),
    ("Dec 2024", "Fructose (mix)", WK16, "*10gL_Glu_Fru*", (14.6, 15.0), 10.0, False),
    ("Dec 2024", "Lactic acid", WK16, "*_LacticAcid*", (18.5, 18.9), 5.0, True),
    ("Dec 2024", "Acetic acid", WK16, "*_Acetate*", (23.5, 23.9), None, True),
    ("Feb 2024", "Glucose", FEB24, "YH_MB_20240225_Protocol1_Glucose.lcb001*", (25.0, 26.0), None, False),
    ("Feb 2024", "Mannitol", FEB24, "YH_MB_20240225_Protocol1_Mannitol.lcb001*", (28.0, 29.0), None, False),
    ("Feb 2024", "Arabinose", FEB24, "YH_MB_20240225_Protocol1_Arabinose.lcb001*", (30.0, 31.0), None, False),
    ("Feb 2024", "Arabitol", FEB24, "YH_MB_20240225_Protocol1_Arabitol.lcb001*", (30.6, 31.6), None, False),
    ("Feb 2024", "Xylitol", FEB24, "YH_MB_20240225_Protocol1_Xylitol.lcb001*", (31.6, 32.6), None, False),
    ("Jan 2025", "Xylitol", JAN25, "Xylitol_Std_1gl*_001.lcd", (16.0, 16.6), 1.0, False),
]

# Published and lab tables. Values are typed from the sources named here; None means the table has no entry.
RT_SOURCES = [
    ("Lab cheat sheet, Jan 2025 (rt_library_protocol1.csv)", "same protocol, single fixed table", {"Glucose": 13.014, "Fructose": 14.306, "Lactic Acid 1": 14.5, "Lactic Acid 2": 16.2, "Xylitol": 16.294, "Succinic acid": 18.1, "Glycerol": 19.616, "Acetic acid": 22.9, "Ethanol": 31.568}),
    ("Banks PhD thesis, Table G-3 (Protocol 1)", "same protocol, standards of Jan-Feb 2025", {"Glucose": 13.0, "Fructose": 14.3, "Xylitol": 16.3, "Glycerol": 19.6, "Acetic acid": 23.0, "Ethanol": 31.5}),
    ("Bio-Rad Bulletin 1928, Table 4.1 (HPX-87H)", "65 C, 5 mM H2SO4, 0.6 mL/min", {"Glucose": 10.16, "Xylose": 10.39, "Fructose": 10.39, "Arabinose": 11.34}),
    ("Phenomenex Rezex RHM H+ chart (HPX-87H equivalent)", "80 C, water, 0.6 mL/min", {"Glucose": 9.62, "Mannitol": 10.59, "Sorbitol": 10.83, "Ribitol": 11.11, "Arabinose": 11.24, "Ribose": 11.42, "Xylitol": 11.78, "Erythritol": 12.14}),
]


def one(folder: Path, pattern: str) -> Path:
    hits = sorted(p for p in folder.glob(pattern) if p.suffix == ".lcd")
    if not hits:
        raise FileNotFoundError(f"{folder}/{pattern}")
    return hits[0]


def load_pda(path: Path):
    """Wavelengths (nm), time (min) and absorbance matrix (uAU) from the PDA 3D stream, or None if there is no PDA data."""
    ole = olefile.OleFileIO(str(path))
    try:
        if not ole.exists("PDA 3D Raw Data/3D Raw Data"):
            return None
        raw = ole.openstream("PDA 3D Raw Data/3D Raw Data").read()
        table = ole.openstream("PDA 3D Raw Data/Wavelength Table").read()
        item = ole.openstream("PDA 3D Raw Data/3D Data Item").read().decode("latin1")
    finally:
        ole.close()
    if not raw:
        return None
    rate = int(re.search(r"<Rate>(\d+)</Rate>", item).group(1))
    n_wl = struct.unpack_from("<i", table, 0)[0]
    wl = np.array(struct.unpack_from(f"<{n_wl}i", table, 4)) / 100.0
    pos, scans = 0, []
    while pos + 24 <= len(raw) and raw[pos:pos + 2] == b"RC":
        _, n, total = struct.unpack_from("<iii", raw, pos + 4)
        p, values, got, end = pos + 24, [], 0, pos + total
        while got < n and p + 2 <= end:
            nb = struct.unpack_from("<h", raw, p)[0]
            p += 2
            if nb <= 0:
                break
            d = _deltas(raw[p:p + nb])
            values.append(d)
            got += len(d)
            p += nb + 2
        scans.append(np.concatenate(values)[:n] if values else np.full(n, np.nan))
        pos = end
    m = np.vstack(scans)
    return wl, np.arange(m.shape[0]) * rate / 60000.0, m


def max_plot_correlation(path: Path) -> float:
    """Decoder check: the maximum over wavelength of each decoded scan against the Max Plot the instrument stored."""
    wl, t, m = load_pda(path)
    ole = olefile.OleFileIO(str(path))
    try:
        raw = ole.openstream("PDA 3D Raw Data/Max Plot").read()
    finally:
        ole.close()
    n = struct.unpack_from("<i", raw, 8)[0]
    pos, parts, got = 24, [], 0
    while got < n and pos + 2 <= len(raw):
        nb = struct.unpack_from("<h", raw, pos)[0]
        pos += 2
        if nb <= 0:
            break
        d = _deltas(raw[pos:pos + nb])
        parts.append(d)
        got += len(d)
        pos += nb + 2
    stored = np.concatenate(parts)[:n]
    k = min(len(stored), m.shape[0])
    return float(np.corrcoef(m.max(axis=1)[:k], stored[:k])[0, 1])


def status(path: Path) -> dict:
    """Median of each status trace of an injection (oven temperature, pump pressure, ...)."""
    ole = olefile.OleFileIO(str(path))
    out = {}
    try:
        for seg in re.findall(r"<DII[^>]*>(.*?)</DII>", _xml(ole, "LSS Raw Data/2D Data Item"), re.S):
            source = _tag(seg, "DSID") or ""
            if ".CH#" in source:
                continue
            stream = f"LSS Raw Data/StatusLog Ch{_tag(seg, 'CN')}"
            if not ole.exists(stream):
                continue
            raw = ole.openstream(stream).read()
            if len(raw) <= 24:
                continue
            m = re.search(r"<CF>([0-9A-Fa-f]{16})</CF>", seg)
            cf = struct.unpack("<d", bytes.fromhex(m.group(1)))[0] if m else 1.0
            n = struct.unpack_from("<i", raw, 8)[0]
            pos, parts, got = 24, [], 0
            while got < n and pos + 2 <= len(raw):
                nb = struct.unpack_from("<h", raw, pos)[0]
                pos += 2
                if nb <= 0:
                    break
                d = _deltas(raw[pos:pos + nb])
                parts.append(d)
                got += len(d)
                pos += nb + 2
            if parts:
                out[str(_tag(seg, "DN"))] = float(np.median(np.concatenate(parts)[:n] * cf))
    finally:
        ole.close()
    return out


def central_area(t: np.ndarray, y: np.ndarray, apex: float) -> float:
    """Area of apex +- CENTRAL above a cubic baseline fitted on the two flanks. t in min, result in signal x s."""
    near = (np.abs(t - apex) >= FLANK[0]) & (np.abs(t - apex) <= FLANK[1])
    mid = np.abs(t - apex) <= CENTRAL
    if near.sum() < 8 or mid.sum() < 3:
        return float("nan")
    k = np.polyfit(t[near] - apex, y[near], 3)
    return float(np.trapezoid(y[mid] - np.polyval(k, t[mid] - apex), t[mid] * 60.0))


def rid_signal(path: Path):
    inj = read_injection(path)
    t, y = chromatogram(path, inj.rid)
    return inj, np.asarray(t, float), np.asarray(y, float)


def fwhm(t: np.ndarray, y: np.ndarray, apex: float) -> float:
    near = (np.abs(t - apex) >= FLANK[0]) & (np.abs(t - apex) <= FLANK[1])
    k = np.polyfit(t[near] - apex, y[near], 3)
    m = np.abs(t - apex) < FLANK[0]
    tt, s = t[m], y[m] - np.polyval(k, t[m] - apex)
    i = j = int(np.argmax(s))
    half = s[i] / 2
    while i > 0 and s[i] > half:
        i -= 1
    while j < len(s) - 1 and s[j] > half:
        j += 1
    left = tt[i] + (half - s[i]) * (tt[i + 1] - tt[i]) / (s[i + 1] - s[i])
    right = tt[j - 1] + (half - s[j - 1]) * (tt[j] - tt[j - 1]) / (s[j] - s[j - 1])
    return float(right - left)


def apex_in(inj, t, y, window) -> tuple[float, float | None]:
    """Apex from the instrument's peak table when it has a peak in the window, otherwise from the raw signal."""
    peaks = [p for p in inj.peaks if window[0] <= p.rt <= window[1]]
    if peaks:
        best = max(peaks, key=lambda p: p.area)
        return best.rt, best.area
    m = (t >= window[0]) & (t <= window[1])
    return float(t[m][np.argmax(y[m])]), None


def uv_ratios(path: Path, rid_area: float, apex: float, delay: float) -> dict:
    pda = load_pda(path)
    if pda is None or not np.isfinite(rid_area) or rid_area <= 0:
        return {f"uv{lam}_ri": np.nan for lam in (195, 200, 210, 220)}
    wl, tp, m = pda
    out = {}
    for lam in (195, 200, 210, 220):
        out[f"uv{lam}_ri"] = central_area(tp, m[:, int(np.argmin(np.abs(wl - lam)))], apex - delay) / rid_area
    return out


def delays() -> pd.DataFrame:
    """RID lags the PDA. Measured on acids, which are strong on both detectors, in each period that is compared by UV."""
    rows = []
    for period, compound, folder, pattern, window, _, acid in SPEC:
        if not acid or period in ("Feb 2024",) or "minor" in compound:
            continue
        path = one(folder, pattern)
        inj, t, y = rid_signal(path)
        pda = load_pda(path)
        if pda is None:
            continue
        wl, tp, m = pda
        apex, _ = apex_in(inj, t, y, window)
        u = m[:, int(np.argmin(np.abs(wl - 210)))]
        look = (tp >= apex - 0.7) & (tp <= apex)
        # Quadratic refinement of the apex on both detectors, because the PDA samples every 0.64 s
        ru = np.abs(t - apex) <= 0.12
        kr = np.polyfit(t[ru] - apex, y[ru], 2)
        apex_r = apex - kr[1] / (2 * kr[0])
        top = tp[look][np.argmax(u[look])]
        us = np.abs(tp - top) <= 0.12
        ku = np.polyfit(tp[us] - top, u[us], 2)
        apex_u = top - ku[1] / (2 * ku[0])
        rows.append({"period": period, "compound": compound, "rid_apex": apex_r, "pda_apex": apex_u, "delay_min": apex_r - apex_u})
    return pd.DataFrame(rows)


def standards(delay_by_period: dict) -> pd.DataFrame:
    rows = []
    for period, compound, folder, pattern, window, level, acid in SPEC:
        path = one(folder, pattern)
        inj, t, y = rid_signal(path)
        apex, table_area = apex_in(inj, t, y, window)
        area = central_area(t, y, apex)
        st = status(path)
        row = {"period": period, "compound": compound, "acid": acid, "level_gl": level, "inj_uL": inj.inj_vol_uL, "rt": apex, "table_area": table_area,
               "central_area": area, "fwhm_min": fwhm(t, y, apex), "oven_C": st.get("Oven Temp."), "pump_bar": st.get("Pump A Pressure"), "file": path.name}
        delay = delay_by_period.get(period)
        row.update(uv_ratios(path, area, apex, delay) if delay else {f"uv{lam}_ri": np.nan for lam in (195, 200, 210, 220)})
        rows.append(row)
    return pd.DataFrame(rows)


def classify(name: str):
    prefixes = {"19012026_Glu_15gpL": "Glc15", "02122025Glu10Fru5": "Glc10Fru5", "04122025Glu10Fru5": "Glc10Fru5", "15122025Fru10Glu5": "Glc5Fru10",
                "08122025Glu5Fru10": "Glc5Fru10 (second run)", "12012026Fru15": "Fru15"}
    for prefix, condition in prefixes.items():
        if name.startswith(prefix) and re.search(r"TP\d+_", name):
            return condition
    return None


def flask_peaks(delay: float) -> pd.DataFrame:
    rows = []
    for path in sorted(FLASK.rglob("*.lcd")):
        condition = classify(path.name)
        if condition is None:
            continue
        inj = read_injection(path)
        glucose = [p for p in inj.peaks if 13.0 <= p.rt <= 13.7 and p.area > 3000]
        glycerol = [p for p in inj.peaks if 19.6 <= p.rt <= 20.0 and p.area > 2000]
        target = [p for p in inj.peaks if 15.95 <= p.rt <= 16.25]
        row = {"condition": condition, "file": path.name, "inj_uL": inj.inj_vol_uL,
               "glucose_rt": max(glucose, key=lambda p: p.area).rt if glucose else np.nan,
               "glycerol_rt": max(glycerol, key=lambda p: p.area).rt if glycerol else np.nan,
               "rt": np.nan, "table_area": 0.0, "central_area": np.nan, "fwhm_min": np.nan}
        if target:
            best = max(target, key=lambda p: p.area)
            _, t, y = rid_signal(path)
            area = central_area(t, y, best.rt)
            row.update(rt=best.rt, table_area=best.area, central_area=area, fwhm_min=fwhm(t, y, best.rt) if area > 0 else np.nan)
            row.update(uv_ratios(path, area, best.rt, delay))
        rows.append(row)
    return pd.DataFrame(rows)


def linear_map(t, g0, y0, g1, y1):
    """Retention time on system 1 of a neutral compound seen at t on system 0, using glucose (g) and glycerol (y) as anchors."""
    return g1 + (t - g0) * (y1 - g1) / (y0 - g0)


def anchored(flask: pd.DataFrame, std: pd.DataFrame, glucose_may: float) -> tuple[pd.DataFrame, dict]:
    may = std[std.period == "May 2025"].set_index("compound").rt
    glycerol_may = float(may["Glycerol (mix L1)"])
    x_feb, x_rezex = arabinose_x(std)
    rows = []
    # Peaks well above noise only: the retention of a very small peak on a sloping baseline is not informative
    for condition, g in flask[(flask.central_area > 2000)].groupby("condition"):
        run = flask[flask.condition == condition]
        g1 = float(run.glucose_rt.median())
        y1 = float(run.glycerol_rt.median())
        row = {"condition": condition, "n": len(g), "glucose_rt": g1, "glycerol_rt": y1, "rt_median": g.rt.median(), "rt_iqr": g.rt.quantile(0.75) - g.rt.quantile(0.25)}
        for compound in ("Arabitol", "Xylitol", "Sorbitol", "Mannitol", "Malic acid"):
            pred = linear_map(float(may[compound]), glucose_may, glycerol_may, g1, y1)
            row[f"pred_{compound}"] = pred
            row[f"delta_{compound}"] = g.rt.median() - pred
        xyl = row["pred_Xylitol"]
        for lab, x in (("Arabinose (x from 25 C)", x_feb), ("Arabinose (x from 80 C)", x_rezex)):
            pred = g1 + x * (xyl - g1)
            row[f"pred_{lab}"] = pred
            row[f"delta_{lab}"] = g.rt.median() - pred
        rows.append(row)
    wk = std[std.period == "Dec 2024"].set_index("compound").rt
    # How far does such a map miss a known compound? Dec 2024 to May 2025 with glucose and mannitol as anchors, arabitol as the test
    predicted = glucose_may + (float(wk["Arabitol"]) - float(wk["Glucose (mix)"])) * (float(may["Mannitol"]) - glucose_may) / (float(wk["Mannitol"]) - float(wk["Glucose (mix)"]))
    check = {"arabitol_may2025_predicted_from_10C_system": predicted, "arabitol_may2025_measured": float(may["Arabitol"]), "miss_min": predicted - float(may["Arabitol"])}
    return pd.DataFrame(rows), check


def arabinose_x(std: pd.DataFrame) -> tuple[float, float]:
    """Position of arabinose between glucose (0) and xylitol (1), at 25 C on this instrument and at 80 C in the Phenomenex chart."""
    feb = std[std.period == "Feb 2024"].set_index("compound").rt
    x_feb = (float(feb["Arabinose"]) - float(feb["Glucose"])) / (float(feb["Xylitol"]) - float(feb["Glucose"]))
    x_rezex = (11.24 - 9.62) / (11.78 - 9.62)
    return x_feb, x_rezex


def arabinose_position(std: pd.DataFrame, glucose_may: float) -> pd.DataFrame:
    feb = std[std.period == "Feb 2024"].set_index("compound").rt
    may = std[std.period == "May 2025"].set_index("compound").rt
    rez = {"Mannitol": 10.59, "Sorbitol": 10.83, "Arabinose": 11.24}
    rows = []
    for c in ("Mannitol", "Sorbitol", "Arabinose", "Arabitol"):
        rows.append({"compound": c,
                     "x_25C_this_instrument": (float(feb[c]) - float(feb["Glucose"])) / (float(feb["Xylitol"]) - float(feb["Glucose"])) if c in feb else np.nan,
                     "x_30C_standards_may2025": (float(may[c]) - glucose_may) / (float(may["Xylitol"]) - glucose_may) if c in may else np.nan,
                     "x_80C_phenomenex": (rez[c] - 9.62) / (11.78 - 9.62) if c in rez else np.nan})
    return pd.DataFrame(rows)


def may_glucose_anchor() -> float:
    """Median retention time of the largest peak between 13.0 and 13.7 min (glucose) in the 1 uL sample injections of
    the May 2025 sequence."""
    rts = []
    for path in sorted(MAY_SEQUENCE.glob("Centrate_*.lcd")):
        inj = read_injection(path)
        if inj.inj_vol_uL != 1.0:
            continue
        p = max((p for p in inj.peaks if 13.0 <= p.rt <= 13.7), key=lambda p: p.area, default=None)
        rts.append(p.rt if p else np.nan)
    return float(pd.Series(rts, dtype=float).median())


def poly_area(t, y, lo, hi, wlo, whi, deg) -> float:
    """Same integration as glycerol_diagnosis.py, so its fructose-corrected areas and these standards share a window."""
    side = ((t >= wlo) & (t < lo)) | ((t > hi) & (t <= whi))
    mid = (t >= lo) & (t <= hi)
    c = np.polyfit(t[side] - lo, y[side], deg)
    return float(np.trapezoid(y[mid] - np.polyval(c, t[mid] - lo), t[mid] * 60))


def arabitol_window_response(folder: Path, pattern: str, level: float) -> float:
    """Area per g/L per uL in the 16.1 min window, shifted with the apex of that standard."""
    inj, t, y = rid_signal(one(folder, pattern))
    apex, _ = apex_in(inj, t, y, (15.8, 16.4))
    shift = apex - 16.10
    return poly_area(t, y, 15.85 + shift, 16.4 + shift, 15.45 + shift, 16.85 + shift, 3) / (level * inj.inj_vol_uL)


def concentration(std: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Approximate only. Sample areas are the fructose-corrected areas of results/glycerol_diagnosis/peak_16p1.csv, and the standards
    are integrated in the same window. The two arabitol standards differ 1.6-fold, so the result is a range."""
    low = arabitol_window_response(P1, "Arabitol_1 g*", 1.0)
    low_01 = arabitol_window_response(P1, "Arabitol_0.1 g*", 0.1)
    high = arabitol_window_response(WK16, "*_Arabitol*", 2.0)
    s = std.set_index(["period", "compound"])
    others = {c: float(s.loc[("May 2025", c), "central_area"]) for c in ("Xylitol", "Sorbitol", "Mannitol")}
    others["Arabitol"] = float(s.loc[("May 2025", "Arabitol"), "central_area"])
    inj, t, y = rid_signal(one(STD26, "20261005_Gly_2gL001*"))
    glycerol_response = poly_area(t, y, 19.45, 20.25, 18.9, 20.8, 2) / (2.0 * inj.inj_vol_uL)
    d = pd.read_csv(RESULTS / "glycerol_diagnosis" / "peak_16p1.csv")
    rows = []
    for condition, g in d.groupby("condition"):
        by_time = g.groupby("time_h").agg(area=("area_16p1", "mean"), glycerol=("glycerol_window", "mean"))
        t_top = by_time.area.idxmax()
        top = by_time.loc[t_top]
        glycerol_gl = top.glycerol / glycerol_response
        rows.append({"condition": condition, "time_h_at_max": t_top, "max_area_per_uL": top.area, "arabitol_gl_low": top.area / high, "arabitol_gl_high": top.area / low,
                     "glycerol_gl_same_time": glycerol_gl,
                     "arabitol_per_glycerol_mol_low": (top.area / high / 152.15) / (glycerol_gl / 92.09) if glycerol_gl > 0.02 else np.nan,
                     "arabitol_per_glycerol_mol_high": (top.area / low / 152.15) / (glycerol_gl / 92.09) if glycerol_gl > 0.02 else np.nan,
                     "last_time_h": by_time.index[-1], "last_area_per_uL": by_time.area.iloc[-1]})
    return pd.DataFrame(rows), {"arabitol_window_area_per_gl_ul_may2025_1gl": low, "arabitol_window_area_per_gl_ul_may2025_0p1gl": low_01, "arabitol_window_area_per_gl_ul_dec2024_2gl": high,
                                "central_area_of_1gl_polyol_standards_may2025": others, "glycerol_window_area_per_gl_ul": glycerol_response}


def model_check() -> pd.DataFrame:
    """Can iJL1143 make these compounds from glucose? Feasibility only: the exit is opened and maximised with growth free."""
    path = CURATED
    if not path.exists():
        return pd.DataFrame()
    import cobra

    m = cobra.io.read_sbml_model(str(path))
    glucose_met = next(x for x in m.metabolites if x.name.lower() in ("d-glucose (extracellular)", "glucose (extracellular)"))
    glucose = next(r for r in glucose_met.reactions if r.boundary)
    carbon_uptakes = [r.id for r in m.boundary if r.lower_bound < 0 and r is not glucose and any("C" in x.elements for x in r.metabolites)]
    rows = []
    for name in ("D-Arabitol (extracellular)", "L-Arabitol (extracellular)", "L-Arabinose (extracellular)", "D-Arabinose (extracellular)", "Xylitol (extracellular)", "D-Mannitol (extracellular)", "Glycerol (extracellular)"):
        met = next((x for x in m.metabolites if x.name == name), None)
        if met is None:
            rows.append({"compound": name, "exit": None, "max_per_glucose_mol": np.nan, "other_carbon_uptakes_closed": len(carbon_uptakes)})
            continue
        exits = [r for r in met.reactions if r.boundary]
        with m:
            for rid in carbon_uptakes:
                m.reactions.get_by_id(rid).lower_bound = 0
            glucose.lower_bound = -1.0
            exit_r = exits[0] if exits else m.add_boundary(met, type="demand")
            exit_r.upper_bound = 1000
            m.objective = exit_r
            value = m.slim_optimize()
        rows.append({"compound": name, "exit": exits[0].id if exits else "none, demand added", "max_per_glucose_mol": value, "other_carbon_uptakes_closed": len(carbon_uptakes)})
    return pd.DataFrame(rows)


def rt_sources() -> pd.DataFrame:
    return pd.DataFrame([{"source": s, "conditions": c, "compound": k, "rt_min": v} for s, c, d in RT_SOURCES for k, v in d.items()])


def figure_data(std: pd.DataFrame, flask: pd.DataFrame, glucose_may: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Traces and UV ratios behind the identity figure, written to CSV so that paper_figures.py draws from results only.
    Standard retention times are mapped onto the Glc 15 run with glucose and glycerol. Each trace is drawn only where its local
    baseline was fitted, and scaled to its own peak."""
    may = std[std.period == "May 2025"].set_index("compound")
    run = flask[flask.condition == "Glc15"]
    g1, y1 = float(run.glucose_rt.median()), float(run.glycerol_rt.median())
    glycerol_may = float(may.loc["Glycerol (mix L1)", "rt"])
    sample = next(FLASK.rglob("19012026_Glu_15gpL_TP8_A3*"))
    traces = [("Glc 15, 48 h", "sample", sample, False),
              ("Arabitol", "polyol", one(P1, "Arabitol_1 g*"), True),
              ("Xylitol", "polyol", one(P1, "Xylitol_1 g*"), True),
              ("Sorbitol", "polyol", one(P1, "Sorbitol_1 g*"), True),
              ("Mannitol", "polyol", one(P1, "Mannitol_1 g*"), True),
              ("Malic acid", "acid", one(P1, "Malic acid_1 g*"), True)]
    windows = {"Glc 15, 48 h": (15.95, 16.25), "Arabitol": (15.9, 16.4), "Xylitol": (16.4, 16.9),
               "Sorbitol": (15.0, 15.3), "Mannitol": (14.6, 15.0), "Malic acid": (14.2, 14.7)}
    rows = []
    for label, kind, path, mapped in traces:
        inj, t, y = rid_signal(path)
        apex, _ = apex_in(inj, t, y, windows[label])
        near = (np.abs(t - apex) >= FLANK[0]) & (np.abs(t - apex) <= FLANK[1])
        k = np.polyfit(t[near] - apex, y[near], 3)
        win = np.abs(t - apex) <= FLANK[1]
        s = y[win] - np.polyval(k, t[win] - apex)
        s = s / s[np.abs(t[win] - apex) < 0.2].max()
        tx = linear_map(t[win], glucose_may, glycerol_may, g1, y1) if mapped else t[win]
        apex_x = float(linear_map(apex, glucose_may, glycerol_may, g1, y1)) if mapped else apex
        rows += [{"trace": label, "kind": kind, "rt_on_flask_run": a, "signal_scaled": b, "apex_on_flask_run": apex_x} for a, b in zip(tx, s)]
    traces_df = pd.DataFrame(rows)
    oct26 = std[std.period == "Oct 2026"].set_index("compound")
    uv = [("Malic acid", "acid", may.loc["Malic acid", "uv210_ri"]), ("Citric acid", "acid", may.loc["Citric acid", "uv210_ri"]),
          ("Acetic acid", "acid", may.loc["Acetic acid", "uv210_ri"]), ("Lactic acid", "acid", oct26.loc["Lactic acid", "uv210_ri"]),
          ("Fructose", "sugar", oct26.loc["Fructose", "uv210_ri"]), ("Glucose", "sugar", oct26.loc["Glucose", "uv210_ri"]),
          ("Glycerol", "polyol", oct26.loc["Glycerol", "uv210_ri"]), ("Mannitol", "polyol", may.loc["Mannitol", "uv210_ri"]),
          ("Sorbitol", "polyol", may.loc["Sorbitol", "uv210_ri"]), ("Xylitol", "polyol", may.loc["Xylitol", "uv210_ri"]),
          ("Arabitol", "polyol", may.loc["Arabitol", "uv210_ri"])]
    uv_rows = [{"label": a, "kind": b, "uv210_ri": float(c), "file": ""} for a, b, c in uv]
    g15 = flask[(flask.condition == "Glc15") & (flask.central_area > 2000)]
    uv_rows += [{"label": "16.1 min peak, Glc 15", "kind": "sample", "uv210_ri": float(v), "file": f} for v, f in zip(g15.uv210_ri, g15.file)]
    return traces_df, pd.DataFrame(uv_rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pd.set_option("display.width", 220, "display.max_columns", 40, "display.max_rows", 200)
    dl = delays()
    delay_by_period = dl.groupby("period").delay_min.mean().to_dict()
    std = standards(delay_by_period)
    glucose_may = may_glucose_anchor()
    flask = flask_peaks(float(np.mean([delay_by_period["May 2025"], delay_by_period["Oct 2026"]])))
    anch, check = anchored(flask, std, glucose_may)
    arabinose = arabinose_position(std, glucose_may)
    conc, responses = concentration(std)
    gem = model_check()
    sources = rt_sources()
    decoder = {p.name: max_plot_correlation(p) for p in (one(P1, "Arabitol_1 g*"), one(STD26, "20261005_LA_2gL001*"), next(FLASK.rglob("19012026_Glu_15gpL_TP8_A3*")))}

    for name, frame in {"delay": dl, "standards": std, "flask_peaks": flask, "anchored_rt": anch, "arabinose_position": arabinose,
                        "concentration": conc, "model_check": gem, "rt_sources": sources}.items():
        frame.to_csv(OUT / f"{name}.csv", index=False)
    traces, uv = figure_data(std, flask, glucose_may)
    traces.to_csv(OUT / "fig_traces.csv", index=False)
    uv.to_csv(OUT / "fig_uv.csv", index=False)

    g15 = flask[(flask.condition == "Glc15") & (flask.central_area > 2000)]
    summary = {
        "decoder_max_plot_r": decoder,
        "delay_min_by_period": delay_by_period,
        "glucose_rt_may2025_anchor": glucose_may,
        "flask_glc15": {"n": len(g15), "rt_median": g15.rt.median(), "uv210_ri_median": g15.uv210_ri.median(), "uv200_ri_median": g15.uv200_ri.median(), "uv195_ri_median": g15.uv195_ri.median(), "fwhm_median": g15.fwhm_min.median()},
        "uv210_expected_if_acid_area_ratio": float(std[(std.period == "Oct 2026") & (std.compound == "Lactic acid")].uv210_ri.iloc[0]),
        "anchored_cross_check": check,
        "responses": responses,
        "xylitol_rt_jan2025_vs_may2025": [float(std[(std.period == "Jan 2025")].rt.iloc[0]), float(std[(std.period == "May 2025") & (std.compound == "Xylitol")].rt.iloc[0])],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=float))
    write_provenance(OUT, params={
        "engine": "none for the identification (raw instrument files); COBRApy feasibility on iJL1143 for model_check only",
        "raw": "FBAVAL_RAW_DIR/lcd", "central_window_min": CENTRAL, "baseline_flanks_min": list(FLANK), "fitted": "nothing; cubic local baselines only",
        "anchors": "glucose and glycerol retention times of each run", "reuses": "results/glycerol_diagnosis/peak_16p1.csv for the re-integrated areas"})

    print(dl.round(3).to_string())
    print(std.drop(columns=["file"]).round(3).to_string())
    print(anch.round(3).T.to_string())
    print(arabinose.round(3).to_string())
    print(conc.round(3).to_string())
    print(gem.round(3).to_string())
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()
