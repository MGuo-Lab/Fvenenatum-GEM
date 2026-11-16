"""Identity of the HPLC products in the batch fermentations, from the raw instrument files.

Reads every RID peak of the fermentation injections and the October 2026 standards. Three questions:
  1. Which raw peak carries the area each workbook calls lactate?
  2. Where does lactic acid elute? The single standard gives the retention times and the response.
  3. Is there lactic acid in the fermentations? The RID signal in the window of its main peak is integrated
     above a straight baseline for every injection, with the mixed standards as the positive control.
Also records per injection the peaks at 14.8 and 16.1 min, which no standard run here explains. The 16.1 min area is
integrated with fbaval.chromatogram, because the peak table is inconsistent on the fructose tail; the table value stays
in peak_16p1_table_area_per_ul.

Raw files: the Shimadzu LabSolutions data files (.lcd) of the instrument, under FBAVAL_RAW_DIR/lcd (see fbaval.paths);
they are not part of this repository. Reader: fbaval.lcd (requires olefile).

Outputs: results/flask_hplc/{peaks.csv, area_match.csv, standards.csv, lactic_control.csv, lactic_window.csv,
         unidentified.csv, overlay.csv, summary.json, provenance.json}

Usage (from analysis/):  python scripts/flask_hplc.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

from fbaval import data
from fbaval.chromatogram import peak_16p1_area
from fbaval.lcd import chromatogram, read_injection
from fbaval.paths import RAW_DIR, RESULTS
from fbaval.provenance import write_provenance

OUT = RESULTS / "flask_hplc"
RAW = RAW_DIR / "lcd" / "TV"
STD = RAW / "HPLC Standards" / "20261005_Standards"
RUNS = {"Glc15": "19012026_Glucose_15gpL", "Glc10Fru5": "02122025GluFru", "Glc5Fru10": "15122025_FruGlu", "Fru15": "15012026_Fru"}
OTHER_RUN = "Glu5Fru10_TP"     # the Glc10Fru5 folder also holds injections of a separate run, with no workbook or online data
LACTIC = (18.92, 19.32)        # window around the main lactic acid peak
FLANK = 0.08                   # min, baseline windows on each side
UNIDENTIFIED = {"peak_14p8": (14.55, 14.95), "peak_16p1": (15.95, 16.25)}


def window_area(t: np.ndarray, y: np.ndarray) -> float:
    left = (t >= LACTIC[0] - FLANK) & (t < LACTIC[0])
    right = (t > LACTIC[1]) & (t <= LACTIC[1] + FLANK)
    mid = (t >= LACTIC[0]) & (t <= LACTIC[1])
    x0, x1, y0, y1 = t[left].mean(), t[right].mean(), y[left].mean(), y[right].mean()
    return float(np.trapezoid(y[mid] - (y0 + (y1 - y0) * (t[mid] - x0) / (x1 - x0)), t[mid] * 60.0))


def signal(path: Path):
    inj = read_injection(path)
    t, y = chromatogram(path, inj.rid)
    return inj, np.asarray(t, float), np.asarray(y, float)


def workbook_areas(cond: str) -> pd.DataFrame:
    wb = openpyxl.load_workbook(data.condition_dir(cond) / "workbook.xlsx", data_only=True)
    rows = list(wb["HPLC"].iter_rows(values_only=True))
    hdr = [str(h or "") for h in rows[0]]
    t_col = hdr.index("Elapsed Time (h)")
    out = []
    for comp in ("Glucose", "Fructose", "Glycerol", "Lactate"):
        cols = [j for j, h in enumerate(hdr) if h.startswith(f"{comp} Peak Area")]
        for r in rows[1:]:
            if isinstance(r[0], str) and r[0].startswith("TP"):
                for fl, j in enumerate(cols, 1):
                    a = r[j] if isinstance(r[j], (int, float)) else 0.0
                    out.append((int(r[0][2:]), float(r[t_col]), fl, comp, float(a)))
    return pd.DataFrame(out, columns=["tp", "time_h", "flask", "compound", "wb_area"])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    peaks, windows, paths = [], [], {}
    for cond, folder in RUNS.items():
        for f in sorted((RAW / "Batch Shake Flask Winter 2025-26" / folder).rglob("*.lcd")):
            m = re.search(r"TP(\d+)_A(\d)", f.name)
            if not m or (cond == "Glc10Fru5" and OTHER_RUN in f.name):
                continue
            paths[f.name] = f
            inj, t, y = signal(f)
            tp, fl = int(m.group(1)), int(m.group(2))
            windows.append({"condition": cond, "tp": tp, "flask": fl, "file": f.name, "inj_uL": inj.inj_vol_uL, "window_area": window_area(t, y)})
            peaks += [{"condition": cond, "tp": tp, "flask": fl, "file": f.name, "inj_uL": inj.inj_vol_uL, "rt": p.rt, "area": p.area, "manual": p.manual} for p in inj.peaks]
    peaks = pd.DataFrame(peaks)
    windows = pd.DataFrame(windows)

    # standards: retention time and response of each single analyte; the mixed series as positive control
    std_rows, control = [], []
    for f in sorted(STD.glob("*.lcd")):
        if not re.search(r"(_LA_2gL|_Gly_2gL|_EtOH_2gL|_Acetic_2gL|EtOH_Gly_LA)", f.name):
            continue
        inj, t, y = signal(f)
        level = float(re.search(r"_(\d+(?:p\d+)?)gL", f.name).group(1).replace("p", "."))
        kind = "mixed" if "EtOH_Gly_LA" in f.name else re.search(r"20261005_([A-Za-z]+)_", f.name).group(1)
        std_rows += [{"standard": kind, "level_gl": level, "file": f.name, "inj_uL": inj.inj_vol_uL, "rt": p.rt, "area": p.area} for p in inj.peaks if p.rt > 12]
        if kind in ("mixed", "LA"):
            control.append({"standard": kind, "level_gl": level, "file": f.name, "inj_uL": inj.inj_vol_uL, "window_area": window_area(t, y),
                            "table_area_19p1": next((p.area for p in inj.peaks if LACTIC[0] < p.rt < LACTIC[1]), np.nan)})
    std = pd.DataFrame(std_rows)
    control = pd.DataFrame(control)
    control["per_gl_ul"] = control.window_area / control.level_gl / control.inj_uL
    # the mixed 1 g/L standard reads low on most analytes (ARIA section 4), so the response is the median over levels
    response = float(control.per_gl_ul.median())
    windows["lactic_upper_gl"] = windows.window_area.clip(lower=0) / response / windows.inj_uL

    # which raw peak carries each workbook area, and the injection the workbook used at each time point
    match_rows, unid = [], []
    for cond in RUNS:
        wb = workbook_areas(cond)
        for (tp, fl), g in wb.groupby(["tp", "flask"]):
            q = peaks[(peaks.condition == cond) & (peaks.tp == tp) & (peaks.flask == fl)]
            if q.empty:
                continue
            refs = g[g.compound.isin(["Glucose", "Fructose", "Glycerol"]) & (g.wb_area > 0)]
            score = {}
            for file, qq in q.groupby("file"):
                score[file] = sum(np.min(np.abs(qq.area - a)) / a for a in refs.wb_area) if len(refs) else 0.0
            best = min(score, key=score.get)
            qb = q[q.file == best]
            for _, r in g[g.wb_area > 0].iterrows():
                k = (qb.area - r.wb_area).abs().idxmin()
                match_rows.append({"condition": cond, "tp": tp, "flask": fl, "compound": r.compound, "wb_area": r.wb_area, "raw_rt": qb.rt[k], "raw_area": qb.area[k],
                                   "relative_difference": abs(qb.area[k] - r.wb_area) / r.wb_area, "manual": bool(qb.manual[k]), "file": best})
            row = {"condition": cond, "tp": tp, "flask": fl, "time_h": float(g.time_h.iloc[0]), "file": best, "inj_uL": float(qb.inj_uL.iloc[0])}
            for name, (a, b) in UNIDENTIFIED.items():
                row[f"{name}_area_per_ul"] = float(qb[qb.rt.between(a, b)].area.sum()) / row["inj_uL"]
            unid.append(row)
    match = pd.DataFrame(match_rows)
    unid = pd.DataFrame(unid).rename(columns={"peak_16p1_area_per_ul": "peak_16p1_table_area_per_ul"})

    # 16.1 min peak: with fructose, the reference is the same flask's first sample after time zero, before any product
    areas = []
    for _, u in unid.iterrows():
        inj, t, y = signal(paths[u.file])
        reference = None
        if data.VINESTOCK[u.condition]["fructose_gl"] > 0:
            first = unid[(unid.condition == u.condition) & (unid.flask == u.flask) & (unid.tp >= 1)].sort_values("tp").iloc[0]
            ref_inj, tr, yr = signal(paths[first.file])
            reference = (tr, yr / ref_inj.inj_vol_uL)
        areas.append(peak_16p1_area(t, y / inj.inj_vol_uL, reference))
    unid["peak_16p1_area_per_ul"] = areas
    blank = unid[(unid.condition == "Fru15") & (unid.tp >= 2)].peak_16p1_area_per_ul
    clean = unid[(unid.condition == "Glc15") & (unid.peak_16p1_table_area_per_ul > 300)]

    # traces for the supplementary figure
    traces = {"Lactic acid standard, 2 g/L": sorted(STD.glob("20261005_LA_2gL001*"))[0],
              "Mixed standard, 0.1 g/L each": sorted(STD.glob("20261005_EtOH_Gly_LA_Ace_0p1gL0_05102026_027*"))[0],
              "Glc 15 fermentation, 48 h": sorted((RAW / "Batch Shake Flask Winter 2025-26" / RUNS["Glc15"]).rglob("*TP8_A2*"))[0]}
    ov = []
    for label, f in traces.items():
        inj, t, y = signal(f)
        m = (t >= 14) & (t <= 21)
        ov.append(pd.DataFrame({"trace": label, "inj_uL": inj.inj_vol_uL, "rt_min": t[m], "signal": y[m] - np.median(y[m])}))
    ov = pd.concat(ov, ignore_index=True)

    la = std[std.standard == "LA"]
    summary = {
        "response_window_area_per_gl_ul": response,
        "lactic_standard_peaks": la.groupby(la.rt.round(1)).area.mean().round(0).to_dict(),
        "lactic_standard_total_area_per_gl_ul": float(la.groupby("file").area.sum().mean() / 2.0 / la.inj_uL.iloc[0]),
        "lactic_upper_gl_max": {c: float(g.lactic_upper_gl.max()) for c, g in windows.groupby("condition")},
        "injections_checked": int(len(windows)),
        "lowest_control": control[control.level_gl == control.level_gl.min()].window_area.mean(),
        "workbook_lactate_raw_rt": {c: [float(g.raw_rt.min()), float(g.raw_rt.max())] for c, g in match[(match.compound == "Lactate") & (match.relative_difference < 0.005)].groupby("condition")},
        "workbook_lactate_exact_matches": int(((match.compound == "Lactate") & (match.relative_difference < 0.005)).sum()),
        "workbook_lactate_values": int((match.compound == "Lactate").sum()),
        "peak_16p1_fru15_blank": {"mean": float(blank.mean()), "sd": float(blank.std()), "min": float(blank.min()), "max": float(blank.max()), "n": int(len(blank))},
        "peak_16p1_area_over_table_glc15_median": float((clean.peak_16p1_area_per_ul / clean.peak_16p1_table_area_per_ul).median()),
        "single_standard_factor_per_gl_ul": {k: float(g.groupby("file").area.max().mean() / 2.0 / g.inj_uL.iloc[0]) for k, g in std[std.standard.isin(["Gly", "EtOH", "Acetic"])].groupby("standard")},
    }
    peaks.to_csv(OUT / "peaks.csv", index=False)
    match.to_csv(OUT / "area_match.csv", index=False)
    std.to_csv(OUT / "standards.csv", index=False)
    control.to_csv(OUT / "lactic_control.csv", index=False)
    windows.to_csv(OUT / "lactic_window.csv", index=False)
    unid.to_csv(OUT / "unidentified.csv", index=False)
    ov.to_csv(OUT / "overlay.csv", index=False)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1))
    write_provenance(OUT, params={"engine": "none (raw instrument files)", "raw": str(RAW.relative_to(RAW_DIR)), "lactic_window_min": LACTIC, "unidentified_windows_min": UNIDENTIFIED,
                                  "fitted_parameters": "none"})
    print(json.dumps(summary, indent=1))
    print(match.groupby(["condition", "compound"]).agg(n=("wb_area", "size"), exact=("relative_difference", lambda s: int((s < 0.005).sum())), rt=("raw_rt", "median")).round(2).to_string())
    print(unid.groupby("condition")[["peak_14p8_area_per_ul", "peak_16p1_table_area_per_ul", "peak_16p1_area_per_ul"]].max().round(0).to_string())


if __name__ == "__main__":
    main()
