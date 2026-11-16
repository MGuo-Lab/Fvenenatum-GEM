"""Figure: initial biomass x0 as the lag parameter of the Curated dFBA in batch fermentation.

Draws figures/figS1_lag.{png,svg,pdf} from results/lag_x0_fit with the shared style of
scripts/paper_figures.py (imported for its constants; the other figures are not redrawn), and adds
the figure to figures/manifest.json.

Usage (from analysis/):  python scripts/lag_x0_figure.py
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("paper_figures", ROOT / "scripts" / "paper_figures.py")
pf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pf)
plt, Line2D, data = pf.plt, pf.Line2D, pf.data

RES = pf.RES / "lag_x0_fit"
NAME = "figS1_lag"
CUR = pf.MODEL_C["curated"]
# one colour per flask, dark (glucose only) to light (fructose only), all from the Curated hue
COND_C = dict(zip(("Glc15", "Glc10Fru5", "Glc5Fru10", "Fru15"), pf.CURATED_RAMP))
SHORT = {"Glc15": "Glc 15", "Glc10Fru5": "Glc 10 + Fru 5", "Glc5Fru10": "Glc 5 + Fru 10", "Fru15": "Fru 15"}
RUN_STY = {"apriori": dict(ls=(0, (1.2, 1.3)), lw=1.0), "shared": dict(ls=(0, (4, 1.6)), lw=1.0), "fit": dict(ls="-", lw=1.3)}


def sci(x: float) -> str:
    """4.6e-4 -> '4.6 × 10⁻⁴' in mathtext."""
    e = int(np.floor(np.log10(x)))
    return f"${x / 10 ** e:.1f}\\times10^{{{e}}}$"


def main() -> None:
    prof = pd.read_csv(RES / "profile.csv")
    fits = pd.read_csv(RES / "fits.csv").set_index("condition")
    ser = pd.read_csv(RES / "series.csv")
    conds = list(data.PRIMARY_VINESTOCK)
    x0_shared = float(fits.x0_shared_gl.iloc[0])
    x0_apr = 0.01

    fig = plt.figure(figsize=(pf.FULL, 118 * pf.MM))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.45, left=0.065, right=0.985, top=0.96, bottom=0.1)
    top = outer[0].subgridspec(1, 3, width_ratios=[1.0, 0.56, 0.9], wspace=0.22)
    bot = outer[1].subgridspec(1, 4, wspace=0.3)

    # (a) timing of the Curated run against x0, with the measured timing and the fitted x0
    ax = fig.add_subplot(top[0, 0]); pf.panel(ax, "a", x=-0.17)
    for c in conds:
        g = prof[(prof.condition == c) & (prof.rule == "glucose_first")]
        ax.plot(g.x0_gl, g.t50_h, color=COND_C[c], lw=1.0, marker="o", ms=2.4, zorder=3)
        g = prof[(prof.condition == c) & (prof.rule == "simultaneous")]
        if len(g):      # mixed flasks with simultaneous uptake (as in Fig. 3): growth at about 2 mu_max
            ax.plot(g.x0_gl, g.t50_h, color=COND_C[c], lw=0.8, ls=(0, (1.2, 1.3)), marker="o", ms=1.8,
                    mfc="white", mew=0.5, alpha=0.9, zorder=2)
        ax.axhline(fits.loc[c, "t50_meas_h"], color=COND_C[c], lw=0.6, ls=(0, (3, 2)), zorder=1)
        ax.plot(fits.loc[c, "x0_fit_gl"], fits.loc[c, "t50_fit_h"], marker="D", ms=4.2, color=COND_C[c],
                mec="white", mew=0.5, ls="none", zorder=5)
    ax.axvline(x0_apr, color=pf.G600, lw=0.7, ls=(0, (1.2, 1.3)))
    ax.axvline(x0_shared, color=pf.MEAS, lw=0.7, ls=(0, (4, 1.6)))
    ax.text(x0_apr, 1.0, " assumed\n 0.01", transform=ax.get_xaxis_transform(), fontsize=6, va="top", color=pf.G700)
    ax.text(x0_shared, 1.0, " shared\n " + sci(x0_shared), transform=ax.get_xaxis_transform(),
            fontsize=6, va="top", color=pf.MEAS)
    ax.set_xscale("log")
    ax.set_xlabel("Initial biomass x₀ (g L$^{-1}$)")
    ax.set_ylabel("Time to half final biomass (h)")
    h = [Line2D([], [], color=COND_C[c], lw=1.0, marker="o", ms=2.4) for c in conds] + \
        [Line2D([], [], color=pf.G700, lw=0.8, ls=(0, (1.2, 1.3)), marker="o", ms=1.8, mfc="white", mew=0.5),
         Line2D([], [], color=pf.G700, lw=0.6, ls=(0, (3, 2))), Line2D([], [], marker="D", ms=4.2, color=pf.G700, ls="none")]
    leg_ax = fig.add_subplot(top[0, 1]); leg_ax.axis("off")
    leg_ax.legend(h, [SHORT[c] for c in conds]
                  + ["mixed, simultaneous\nuptake", "measured\n(mean of 3 replicates)", "fitted x₀"],
                  loc="center left", bbox_to_anchor=(-0.12, 0.5), fontsize=6, handlelength=2.0, labelspacing=0.5)

    # (b) timing error per condition: assumed x0, one shared x0, x0 per condition
    ax = fig.add_subplot(top[0, 2]); pf.panel(ax, "b", x=-0.2)
    w = 0.26
    xs = np.arange(len(conds))
    for k, (col, lab, sty) in enumerate([("timing_err_apriori_h", "assumed x₀ (0.01 g L$^{-1}$)", dict(color="white", ec=CUR, hatch="////")),
                                         ("timing_err_shared_h", f"one shared x₀", dict(color=pf.CURATED_RAMP[3], ec=CUR)),
                                         ("timing_err_fit_h", "x₀ per condition", dict(color=CUR, ec=CUR))]):
        ax.bar(xs + (k - 1) * w, fits.loc[conds, col], width=w, lw=0.6, label=lab, **sty)
    ax.axhline(0, color=pf.G800, lw=0.6)
    ax.set_xticks(xs, [SHORT[c].replace(" + ", "\n+ ") for c in conds])
    ax.set_ylim(-37, 28)
    ax.set_xlabel("Condition (sugars, g L$^{-1}$)")
    ax.set_ylabel("Timing error (h)")
    ax.legend(loc="upper right", fontsize=6, handlelength=1.6, borderaxespad=0.3)

    # (c–f) biomass in each flask: measured, and the Curated run at each x0, with no time shift
    for j, c in enumerate(conds):
        ax = fig.add_subplot(bot[0, j]); pf.panel(ax, "cdef"[j], x=-0.3)
        pres = data.load_vinestock_presens(c).dropna(subset=["biomass_gl"])
        fin = data.load_vinestock_final(c).biomass_gl
        for _, d in pres.groupby("flask"):
            ax.plot(d.time_h, d.biomass_gl, color=pf.G600, lw=0.45, zorder=2)
        t_end = pres.time_h.max()
        ax.errorbar(t_end, fin.mean(), yerr=fin.std(ddof=1), fmt="o", color=pf.MEAS, ms=3.2, mec="white", mew=0.4,
                    elinewidth=0.7, capsize=0, zorder=6)
        for run in ("apriori", "shared", "fit"):
            s = ser[(ser.condition == c) & (ser.run == run)]
            ax.plot(s.time_h, s.biomass_gl, color=CUR, zorder=4, **RUN_STY[run])
        ax.set_xlim(0, np.ceil((t_end + 6) / 10) * 10); ax.set_ylim(0, 10.2)
        ax.set_title(pf.FLASK[c], pad=3)
        ax.set_xlabel("Time (h)")
        if j == 0:
            ax.set_ylabel("Biomass (g L$^{-1}$)")
        ax.text(0.03, 0.97, "x₀ = " + sci(fits.loc[c, "x0_fit_gl"]) + " g L$^{-1}$", transform=ax.transAxes,
                fontsize=6, va="top", color=pf.INK)
    h = [Line2D([], [], color=CUR, **RUN_STY[r]) for r in ("apriori", "shared", "fit")] + \
        [Line2D([], [], color=pf.G600, lw=0.6), Line2D([], [], marker="o", ls="", color=pf.MEAS, ms=3.2)]
    fig.legend(h, ["Curated, assumed x₀ (0.01 g L$^{-1}$)", "Curated, one shared x₀", "Curated, x₀ per condition",
                   "Biomass, online (3 replicates)", "Final biomass, dry weight"],
               loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.035), fontsize=6.0, handlelength=2.0, columnspacing=1.0)

    pf.save(fig, NAME, ["lag_x0_fit"])
    man_p = pf.OUT / "manifest.json"
    man = json.loads(man_p.read_text()) if man_p.exists() else {}
    man[NAME] = pf.MANIFEST[NAME]
    man_p.write_text(json.dumps(man, indent=1))
    print("wrote", NAME, pf.MANIFEST[NAME])


if __name__ == "__main__":
    main()
