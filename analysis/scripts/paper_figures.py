"""Manuscript figures and tables, drawn only from the result files in results/ and the model files in model/.

One style for every figure: Liberation Sans (metric-compatible with Arial),
7 pt text, 6.5 pt ticks, 8 pt bold panel letters, 0.8 pt lines, Okabe-Ito
colours with a fixed meaning (models, measured data, compounds, pH-model
versions). Figures are written as PNG (600 dpi), SVG and PDF to figures/;
tables go to figures/tables/ as CSV. figures/manifest.json records which
result folders (and their provenance commits) each figure reads.
Supplementary Fig. S1 is drawn by lag_x0_figure.py. Fig. 1b needs the template
model, downloaded by gem/download_data.py.

FBAVAL_FIGURES redirects the output folder.

Usage (from analysis/):  python scripts/paper_figures.py
"""

from __future__ import annotations

import csv
import json
import os
import re
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from matplotlib.patches import FancyBboxPatch, Patch

from fbaval import data, names, simulate
from fbaval.paths import CURATED, CURATION_DIR, DRAFT, FIGURES, GEM_INPUTS, RESULTS, TEMPLATE

RES = RESULTS
OUT = Path(os.environ.get("FBAVAL_FIGURES", str(FIGURES)))
TAB = OUT / "tables"

MM = 1 / 25.4
FULL, HALF = 180 * MM, 88 * MM
QO2 = r"$q_{\mathrm{O_2,max}}$"


def tex(label: str) -> str:
    """Display names (fbaval.names) with matplotlib subscripts."""
    return label.replace("q_O₂,max", QO2).replace("NH₄⁺", r"NH$_4^+$").replace("O₂", r"O$_2$")


# One palette for every figure: one blue and one orange (Okabe & Ito) with their tints, grey and black, no other hue.
# In model comparisons blue is Curated and orange is Draft; each figure's legend states what the colours mean.
INK, G800, G700, G600, G500, G400, G300, G200, G100 = ("#1A1A1A", "#333333", "#4D4D4D", "#6E6E6E", "#9E9E9E",
                                                       "#BDBDBD", "#D9D9D9", "#EDEDED", "#F5F5F5")
MODEL_C = {"draft": "#D55E00", "curated": "#0072B2"}            # reserved for the two models
MODEL_TINT = {"draft": "#E8B088", "curated": "#88B8D8"}         # 45% tints of the model colours
CURATED_RAMP = ("#003F62", "#0072B2", "#56A3D0", "#BFDCEC")     # shades of the Curated blue, for variants of that model
MODELS = ["draft", "curated"]                                   # the paper compares these two
MODEL_LAB = {k: tex(v[0]) for k, v in names.MODEL.items()}
MODEL_SHORT = {k: tex(v[1]) for k, v in names.MODEL.items()}
MODEL_TICK = {"draft": "Draft", "curated": "Curated"}
MEAS = INK                                                      # measured mean or end point; replicates use G600
HATCH = {"hatch": "////"}
CMP = {"glycerol": MODEL_C["curated"], "peak16": MODEL_C["draft"], "biomass": MODEL_C["curated"], "products": MODEL_TINT["curated"],
       "ethanol": MODEL_C["draft"], "byproducts": "hatch", "co2_ferm": MODEL_TINT["draft"], "co2": G200}
PHV = {"v0": MODEL_C["draft"], "v1": MODEL_TINT["curated"], "v2": CURATED_RAMP[2], "v3": MODEL_C["curated"]}   # orange start, blue as effects are added
# outcome categories (agreement, forced/possible/none); never used for models
STATUS = {"good": G700, "partial": G400, "neutral": G200, "bad": "hatch", "none": G100}


def fill(col: str) -> dict:
    """Bar or patch style: solid grey, or a hatch where no further grey can be told apart."""
    if col in HATCH:
        return dict(facecolor="white", edgecolor=MODEL_C["draft"], hatch=HATCH[col], lw=0)
    return dict(color=col, lw=0)


SUGAR_MARK = {"glucose": "o", "fructose": "s", "sucrose": "D", "xylose": "^"}
FLASK = {c: names.flask(c) for c in names.FLASK_COND}          # 'Glc 15'
PH_LAB = {k: tex(v[1]) for k, v in names.PH.items()}

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "legend.frameon": False,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "lines.linewidth": 0.9, "lines.markersize": 3.2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": G800,
    "axes.labelcolor": INK, "xtick.color": G800, "ytick.color": G800,
    "hatch.linewidth": 0.5, "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.dpi": 600, "figure.dpi": 150,
    "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic",
    "mathtext.bf": "Liberation Sans:bold", "mathtext.cal": "Liberation Sans", "mathtext.sf": "Liberation Sans",
})
MANIFEST: dict[str, dict] = {}
MANIFEST_NOTES: dict[str, dict] = {}


def panel(ax, letter, x=-0.16, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


def save(fig, name, sources):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    prov = {}
    for s in sources:
        pj = RES / s / "provenance.json"
        if pj.exists():
            r = json.loads(pj.read_text())["repos"]
            prov[s] = {k: (v["commit"][:7], v["dirty"]) for k, v in r.items()}
    MANIFEST[name] = {"sources": sources, "provenance": prov, **MANIFEST_NOTES.get(name, {})}


def jload(p):
    return json.loads((RES / p).read_text())


def model_counts(model_file):
    """Metabolites and genes of a model file."""
    m = simulate.load(str(model_file))
    return {"metabolites": len(m.metabolites), "genes": len(m.genes)}


# ── Figure 1: reconstruction, orthology, GPR evidence ─────────────────────────

def orthology_stats():
    """Template genes of the CBS 513.88 iJB1325 and their ortholog support."""
    import cobra
    genes = {g.id for g in cobra.io.read_sbml_model(str(TEMPLATE)).genes}
    rows = [r for r in csv.DictReader(open(GEM_INPUTS / "orthologs_cbs513.tsv"), delimiter="\t") if r["template_gene"] in genes]
    c = Counter(r["method"] for r in rows)
    ident = [float(r["pident_fwd"]) for r in rows if r["method"].startswith("RBH") and r["pident_fwd"]]
    return {"n": len(genes), "mapped": len(rows), "methods": c, "ident": ident}


def gpr_classes(model_file):
    """Reactions with an EC number and a GPR: does the gene's UniProt EC match the reaction EC,
    share its first three levels, differ, or is the gene unannotated?"""
    import curate_model as cm
    ecs = {}
    with open(GEM_INPUTS / "fven_uniprot_ec.tsv", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            ecs[r["Entry"]] = {e.strip() for e in (r["EC number"] or "").split(";") if re.fullmatch(r"\d+\.\d+\.\d+\.\d+", e.strip())}
    m = simulate.load(str(model_file))
    out = Counter()
    for r in m.reactions:
        rec = set(cm.reaction_ecs(r))
        if not rec or not r.genes:
            continue
        gec = set().union(*(ecs.get(g.id, set()) for g in r.genes))
        if not gec:
            out["gene without EC"] += 1
        elif gec & rec:
            out["EC matches"] += 1
        elif any(a.rsplit(".", 1)[0] == b.rsplit(".", 1)[0] for a in rec for b in gec):
            out["related EC (same x.y.z)"] += 1
        else:
            out["different EC"] += 1
    return out


def stage_counts():
    log = list(csv.DictReader(open(CURATION_DIR / "curation_log.csv")))
    c = Counter((r["stage"], r["action"]) for r in log)
    return {"draft_removed": c[("1_draft", "removed")], "ec_rescue": c[("2_ec_rescue", "added")],
            "gapfill": c[("2_gapfill", "added")] + c[("4_task_recheck", "added")],
            "exports_closed": c[("3_secretion", "export_closed")],
            "thermo_v1": sum(v for (s, a), v in c.items() if s == "4_thermo" and a.startswith("direction")),
            "redox": sum(v for (s, a), v in c.items() if s == "4_redox"),
            "thermo_v2": sum(v for (s, a), v in c.items() if s == "6_v2_thermo" and a.startswith("direction")),
            "acid": c[("6_v2_acid_export", "export_charged")]}


def fig1(qc):
    st, orth = stage_counts(), orthology_stats()
    fig = plt.figure(figsize=(FULL, 98 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.95, 1], width_ratios=[1, 1], hspace=0.42, wspace=0.42, left=0.07, right=0.99, top=0.97, bottom=0.1)
    ax = fig.add_subplot(gs[0, :]); ax.set_axis_off(); ax.set_xlim(-1.2, 101.2); ax.set_ylim(0, 30)
    panel(ax, "a", x=-0.02, y=0.97)
    W, Y0, H, GAP = 17.0, 11.5, 18.2, 3.75
    PAD, CLEAR = 0.2, 0.4                                       # box padding; space between an arrow and a box border
    boxes = [
        (0.0, "Template", "iJB1325 (A. niger)\nCBS 513.88 version\n2,320 reactions\n1,373 genes"),
        (W + GAP, "Orthology", "RBH blastp\n+ OrthoFinder\n(10 fungi)\n" + f"{orth['mapped']:,} of {orth['n']:,}\ngenes mapped"),
        (2 * (W + GAP), "Stages 1–5", f"orthology −{st['draft_removed']}, EC +{st['ec_rescue']}\ngap-filling +{st['gapfill']}\n"
                                    f"{st['exports_closed']} exchanges closed\nredox, directionality\nbiomass, NGAM"),
        (3 * (W + GAP), "Stage 6", f"directionality of all\nreactions ({st['thermo_v2']})\norganic-acid export\ncost ({st['acid']} transports)\nGTP/GDP carrier"),
        (4 * (W + GAP), "Curated", f"iJL1143\n{qc['curated']['reactions']:,} reactions\n{qc['curated']['genes']:,} genes\n{qc['curated']['open_exports']} open exchanges"),
    ]
    for x, title, body in boxes:
        col = MODEL_C["curated"] if title == "Curated" else G700
        ax.add_patch(FancyBboxPatch((x, Y0), W, H, boxstyle="round,pad=0.2,rounding_size=1.0", fc="white", ec=col, lw=0.9))
        ax.text(x + W / 2, Y0 + H - 2.4, title, ha="center", va="center", fontweight="bold", color=INK, fontsize=6.8)
        ax.text(x + W / 2, Y0 + H - 4.5, body, ha="center", va="top", fontsize=5.9, color=INK, linespacing=1.3)
    for k in range(1, 5):
        left = k * (W + GAP)                                    # left edge of the box the arrow points to
        ax.annotate("", xy=(left - PAD - CLEAR, Y0 + H / 2), xytext=(left - GAP + PAD + CLEAR, Y0 + H / 2),
                    arrowprops=dict(arrowstyle="-|>", lw=0.8, color=G700, shrinkA=0, shrinkB=0, mutation_scale=7))
    LX = 2 * (W + GAP)
    DW = 2 * W + GAP
    ax.add_patch(FancyBboxPatch((LX, 0.3), DW, 8.6, boxstyle="round,pad=0.2,rounding_size=1.0",
                                fc="white", ec=MODEL_C["draft"], lw=0.9, ls=(0, (3, 2))))
    ax.text(LX + DW / 2, 6.9, "Draft", ha="center", va="center", fontweight="bold", color=INK, fontsize=6.8)
    ax.text(LX + DW / 2, 3.0, f"all template reactions; biomass, NGAM as in Curated\n"
            f"{qc['draft']['reactions']:,} reactions · {qc['draft']['genes']:,} genes · {qc['draft']['open_exports']} open exchanges",
            ha="center", va="center", fontsize=5.9, linespacing=1.3)
    ax.annotate("", xy=(LX - PAD - CLEAR, 4.7), xytext=(W + GAP + W / 2, Y0 - PAD - CLEAR),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color=MODEL_C["draft"], mutation_scale=7, shrinkA=0, shrinkB=0,
                                connectionstyle="angle,angleA=-90,angleB=180,rad=0"))

    ax = fig.add_subplot(gs[1, 0]); panel(ax, "b")
    ax.hist(orth["ident"], bins=np.arange(40, 102.5, 2.5), color=MODEL_C["curated"], edgecolor="white", lw=0.4)
    ax.set_xlabel("Sequence identity of reciprocal best hits (%)"); ax.set_ylabel("Template genes")
    ax.set_xlim(40, 100); ax.set_ylim(0, ax.get_ylim()[1] * 1.55)
    mc = orth["methods"]
    txt = (f"{orth['mapped']:,} of {orth['n']:,} template genes mapped ({100 * orth['mapped'] / orth['n']:.1f}%)\n"
           f"RBH and OrthoFinder agree: {mc.get('RBH+OF', 0)}\nOrthoFinder only: {mc.get('OF', 0)}\n"
           f"RBH only: {mc.get('RBH', 0)} · the two disagree: {mc.get('RBH(OF-disagree)', 0)}")
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", ha="left", fontsize=6, linespacing=1.35)

    ax = fig.add_subplot(gs[1, 1]); panel(ax, "c")
    # each row in the colour of its model: full colour, tint, hatch, then neutral for "no EC"
    def seg(m, level):
        if level == "hatch":
            return dict(facecolor="white", edgecolor=MODEL_C[m], hatch=HATCH["hatch"], lw=0)
        return dict(color={"full": MODEL_C[m], "tint": MODEL_TINT[m], "none": G200}[level], edgecolor="white", lw=0.4)
    cats = [("EC matches", "full"), ("related EC (same x.y.z)", "tint"), ("different EC", "hatch"), ("gene without EC", "none")]
    G = {"draft": gpr_classes(DRAFT), "curated": gpr_classes(CURATED)}
    for i, mk in enumerate(("draft", "curated")):
        left = 0
        for cat, col in cats:
            v = G[mk][cat]
            ax.barh(i, v, left=left, height=0.55, **seg(mk, col))
            if v > 35:
                ax.text(left + v / 2, i, str(v), ha="center", va="center", fontsize=6,
                        color="white" if col == "full" else INK,
                        bbox=dict(facecolor="white", edgecolor="none", pad=0.6) if col == "hatch" else None)
            left += v
    ax.set_yticks(range(2), [MODEL_SHORT["draft"], MODEL_SHORT["curated"]]); ax.invert_yaxis()
    ax.set_xlabel("Reactions with an EC number and a gene association")
    ax.legend([(Patch(**seg("draft", c)), Patch(**seg("curated", c))) for _, c in cats], [l for l, _ in cats],
              handler_map={tuple: HandlerTuple(ndivide=None, pad=0.2)}, handlelength=2.6, loc="upper center",
              bbox_to_anchor=(0.5, -0.32), ncol=2, columnspacing=1)
    save(fig, "fig1_reconstruction", ["model_versions_qc"])
    return st, orth, G


# ── Figure 2: what curation changes ───────────────────────────────────────────

def fig2(qc):
    M = MODELS
    fig, axs = plt.subplots(1, 2, figsize=(84 * MM, 54 * MM), gridspec_kw={"wspace": 0.55, "width_ratios": [1, 1], "left": 0.14, "right": 0.98, "top": 0.89, "bottom": 0.16})
    fig.legend([Patch(color=MODEL_C[m]) for m in M], [MODEL_SHORT[m] for m in M], loc="upper center",
               bbox_to_anchor=(0.5, 1.1), ncol=2, handlelength=1, columnspacing=1.4)
    ax = axs[0]; panel(ax, "a", x=-0.42)
    vals = [qc[m]["max_atp_per_glucose"] for m in M]
    ax.bar(range(2), vals, color=[MODEL_C[m] for m in M], width=0.6)
    for i, (m, v) in enumerate(zip(M, vals)):
        ax.text(i, v - 1.2, f"{v:.1f}", ha="center", va="top", fontsize=6.5, color="white")
        if qc[m]["energy_leaks"]:
            ax.text(i, v + 1.2, f"{len(qc[m]['energy_leaks'])} cycle\ntests fail", ha="center", va="bottom", fontsize=5.8,
                    color=INK, linespacing=1.1)
    ax.axhline(32, color=G600, lw=0.7, ls=(0, (3, 2)), zorder=0)
    ax.text(1.0, 33.0, "theoretical max.\n32 (P/O 2.5)", fontsize=5.6, ha="center", va="bottom", color=G600, linespacing=1.1)
    ax.set_xticks(range(2), [MODEL_TICK[m] for m in M], fontsize=6.5); ax.set_ylim(0, 60); ax.set_xlim(-0.6, 1.6)
    ax.set_ylabel("Maximum ATP yield (mol/mol glucose)")

    ax = axs[1]; panel(ax, "b", x=-0.42)
    S = ["glucose", "fructose"]
    w = 0.36
    for k, m in enumerate(M):
        ax.bar(np.arange(len(S)) + (k - 0.5) * w, [qc[m]["yield_gg"][s] for s in S], width=w, color=MODEL_C[m])
    rep = pd.read_csv(RES / "flask_products" / "replicates.csv")
    meas = {"glucose": rep[rep.condition == "Glc15"].yield_gg.mean(), "fructose": rep[rep.condition == "Fru15"].yield_gg.mean()}
    for s_, y in meas.items():
        ax.plot([S.index(s_) - 0.42, S.index(s_) + 0.42], [y, y], color=MEAS, lw=0.9, ls=(0, (2, 1.4)))
    ax.legend([Line2D([], [], color=MEAS, lw=0.9, ls=(0, (2, 1.4)))], ["measured yield"], loc="upper right")
    ax.set_xticks(range(len(S)), [s_.capitalize() for s_ in S]); ax.set_ylim(0, 0.8)
    ax.set_ylabel("Aerobic biomass yield (g g$^{-1}$)")
    save(fig, "fig2_curation_effects", ["model_versions_qc", "flask_products"])


# ── Figure 3: batch fermentation time courses, final biomass, stepwise error ──

STEPS = [("base", "Base"), ("sugar_order", "+ sequential\nuptake"), ("lag_shared", "+ one\nlag time"), ("lag_each", "+ lag time\nper condition")]


def fig3():
    met = pd.read_csv(RES / "flask_dynamics" / "metrics.csv")
    ser = pd.read_csv(RES / "flask_dynamics" / "series.csv")
    summ = pd.read_csv(RES / "flask_dynamics" / "summary.csv")
    M = MODELS
    STY = {"draft": dict(lw=1.1, zorder=3), "curated": dict(lw=1.3, zorder=4)}
    SLS = {"glucose": "-", "fructose": (0, (1.0, 1.1))}
    SMK = {"glucose": dict(marker="o", mfc=MEAS), "fructose": dict(marker="s", mfc="white")}
    fig = plt.figure(figsize=(FULL, 150 * MM))
    # Margins kept small so that the figure prints at its drawn size in the full text width.
    outer = fig.add_gridspec(2, 1, height_ratios=[1.9, 0.95], hspace=0.3, left=0.075, right=0.985, top=0.975, bottom=0.075)
    top = outer[0].subgridspec(2, 4, height_ratios=[1, 0.85], hspace=0.1, wspace=0.32)
    bot = outer[1].subgridspec(1, 3, width_ratios=[1.1, 1.1, 1.0], wspace=0.42)
    for j, c in enumerate(data.PRIMARY_VINESTOCK):
        ax = fig.add_subplot(top[0, j]); panel(ax, "abcd"[j], x=-0.26)
        axs = fig.add_subplot(top[1, j], sharex=ax)
        pres = data.load_vinestock_presens(c).dropna(subset=["biomass_gl"])
        fin = data.load_vinestock_final(c).biomass_gl
        hplc = data.load_vinestock_hplc(c)
        for _, d in pres.groupby("flask"):
            ax.plot(d.time_h, d.biomass_gl, color=G600, lw=0.45, zorder=2)
        t_end = pres.time_h.max()
        ax.errorbar(t_end, fin.mean(), yerr=fin.std(ddof=1), fmt="o", color=MEAS, ms=3.4, mec="white", mew=0.4,
                    elinewidth=0.7, capsize=0, zorder=6)
        sugars = [sg for sg in ("glucose", "fructose") if data.VINESTOCK[c][f"{sg}_gl"] > 0]
        for m in M:
            sm = ser[(ser.model == m) & (ser.condition == c) & (ser.setting == "lag_each")]
            ax.plot(sm.time_h, sm.biomass_gl, color=MODEL_C[m], **STY[m])
            for sg in sugars:
                axs.plot(sm.time_h, sm[f"{sg}_gl"], color=MODEL_C[m], ls=SLS[sg], **STY[m])
        for sg in sugars:
            h = hplc[hplc.compound == sg].groupby("time_h").conc_gl.agg(["mean", "min", "max"]).reset_index()
            axs.errorbar(h.time_h, h["mean"], yerr=[h["mean"] - h["min"], h["max"] - h["mean"]], ls="none", color=MEAS,
                         ms=2.8, mec=MEAS, mew=0.6, elinewidth=0.5, capsize=0, zorder=6, **SMK[sg])
        ax.set_xlim(0, np.ceil((t_end + 6) / 10) * 10); ax.set_ylim(0, 10.2)
        axs.set_ylim(0, 19)
        ax.tick_params(labelbottom=False)
        ax.set_title(FLASK[c], pad=3)
        axs.set_xlabel("Time (h)")
        if j == 0:
            ax.set_ylabel("Biomass (g L$^{-1}$)")
            axs.set_ylabel("Sugar (g L$^{-1}$)")
            ax.legend([Line2D([], [], color=MODEL_C[m], lw=STY[m]["lw"]) for m in M] +
                      [Line2D([], [], color=G600, lw=0.6), Line2D([], [], marker="o", ls="", color=MEAS, ms=3.4)],
                      [MODEL_SHORT[m] for m in M] + ["Online biomass\n(3 replicates)", "Final dry weight"],
                      loc="upper left", bbox_to_anchor=(0.0, 1.02), borderaxespad=0, handlelength=1.4, handletextpad=0.5,
                      labelspacing=0.3, fontsize=5.8)
        if j == 1:
            axs.legend([Line2D([], [], color=G700, ls=SLS["glucose"], lw=1.0, marker="o", mfc=MEAS, mec=MEAS, ms=2.8),
                        Line2D([], [], color=G700, ls=SLS["fructose"], lw=1.0, marker="s", mfc="white", mec=MEAS, ms=2.8)],
                       ["Glucose: HPLC, model", "Fructose: HPLC, model"],
                       loc="upper left", bbox_to_anchor=(0.02, 1.02), borderaxespad=0, handlelength=1.8, handletextpad=0.5,
                       labelspacing=0.3, fontsize=5.8)

    # e, f: sequential uptake with nothing fitted, the setting the text reports for final biomass.
    fb = met[met.setting == "sugar_order"].set_index(["model", "condition"])
    conds = list(data.PRIMARY_VINESTOCK)

    def discrepancy_bars(ax, vals, fmt, ylabel, ymax):
        w = 0.38
        xpos = np.r_[np.arange(4), 4.35]
        for k, m in enumerate(M):
            v = list(vals[m]) + [float(np.mean(vals[m]))]
            bars = ax.bar(xpos + (k - 0.5) * w, v, width=w, color=MODEL_C[m], label=MODEL_SHORT[m])
            for b_, y in zip(bars, v):
                # The lower Curated bar stands beside the taller Draft bar, so its value starts at the bar's left edge and runs into the gap.
                x_ = (b_.get_x() + b_.get_width() / 2, "center") if m == "draft" else (b_.get_x() + 0.02, "left")
                ax.text(x_[0], y + 0.015 * ymax, fmt.format(y), ha=x_[1], va="bottom", fontsize=5.6, color=G800)
        ax.axvline(3.67, color=G400, lw=0.6, ls=(0, (2, 2)))
        ax.set_xticks(xpos, ["Glc 15", "Glc 10\n+ Fru 5", "Glc 5\n+ Fru 10", "Fru 15", "Mean"], fontsize=6.0)
        ax.set_ylabel(ylabel); ax.set_ylim(0, ymax); ax.set_xlim(-0.55, 4.85)
        ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.04), borderaxespad=0, handlelength=1.2, labelspacing=0.3, fontsize=5.8)

    ax = fig.add_subplot(bot[0, 0]); panel(ax, "e", x=-0.17)
    discrepancy_bars(ax, {m: [abs(fb.err_final_pct[(m, c)]) for c in conds] for m in M}, "{:.1f}",
                     "Absolute discrepancy\nin final biomass (%)", 25)
    ax = fig.add_subplot(bot[0, 1]); panel(ax, "f", x=-0.17)
    discrepancy_bars(ax, {m: [abs(fb.x_final_gl[(m, c)] - fb.x_final_meas_gl[(m, c)]) for c in conds] for m in M}, "{:.2f}",
                     "Absolute discrepancy\nin final biomass (g L$^{-1}$)", 1.85)
    ax.set_yticks([0, 0.5, 1.0, 1.5])

    ax = fig.add_subplot(bot[0, 2]); panel(ax, "g", x=-0.2)
    for m in M:
        q = summ[summ.model == m].set_index("setting")
        ax.plot(range(4), [q.rmse_x[k] for k, _ in STEPS], color=MODEL_C[m], marker="o", ms=3, lw=1.0)
        ax.plot(range(4), [q.rmse_s[k] for k, _ in STEPS], color=MODEL_C[m], marker="s", ms=2.8, lw=0.9, ls=(0, (3, 1.6)), mfc="white")
    ax.set_xticks(range(4), [l for _, l in STEPS], fontsize=6.0); ax.set_xlim(-0.3, 3.3); ax.set_ylim(0, None)
    ax.set_ylabel("Time-course RMSE (g L$^{-1}$)")
    ax.legend([Line2D([], [], color=MODEL_C[m], lw=1.1) for m in M] +
              [Line2D([], [], color=G700, marker="o", ms=3, lw=1.0),
               Line2D([], [], color=G700, marker="s", ms=2.8, lw=0.9, ls=(0, (3, 1.6)), mfc="white")],
              [MODEL_SHORT[m] for m in M] + ["Biomass", "Sugars"],
              loc="upper right", bbox_to_anchor=(1.0, 1.04), borderaxespad=0, handlelength=2.2, labelspacing=0.3, fontsize=5.8)
    save(fig, "fig3_fermentation_time_courses", ["flask_dynamics"])


# ── Figure 4: transient products and dissolved oxygen ─────────────────────────

def fig_products():
    prod = pd.read_csv(RES / "flask_products" / "products.csv")
    on = pd.read_csv(RES / "flask_products" / "online.csv")
    rep = pd.read_csv(RES / "flask_products" / "replicates.csv")
    unid = pd.read_csv(RES / "flask_hplc" / "unidentified.csv")
    unid["area"] = unid.peak_16p1_area_per_ul.clip(lower=0)      # the fructose blank scatters around zero
    C16 = CMP["peak16"]
    fig = plt.figure(figsize=(FULL, 118 * MM))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 0.8, 0.9], hspace=0.13, wspace=0.27, left=0.09, right=0.985, top=0.97, bottom=0.17)
    top = unid.area.max() * 1.08
    for j, c in enumerate(data.PRIMARY_VINESTOCK):
        ax = fig.add_subplot(gs[0, j]); panel(ax, "abcd"[j], x=-0.24)
        axp = fig.add_subplot(gs[1, j], sharex=ax)
        axo = fig.add_subplot(gs[2, j], sharex=ax)
        q = prod[(prod.condition == c) & (prod.compound == "glycerol")].groupby("time_h").conc_gl.agg(["mean", "min", "max"]).reset_index()
        ax.errorbar(q.time_h, q["mean"], yerr=[q["mean"] - q["min"], q["max"] - q["mean"]], color=CMP["glycerol"], marker="o", ms=2.8, lw=0.9, elinewidth=0.5, capsize=0)
        u = unid[unid.condition == c].groupby("time_h").area.agg(["mean", "min", "max"]).reset_index()
        axp.errorbar(u.time_h, u["mean"], yerr=[u["mean"] - u["min"], u["max"] - u["mean"]], color=C16, marker="s", ms=2.6, mfc="white", mew=0.6, lw=0.9, elinewidth=0.5, capsize=0)
        for _, d in on[on.condition == c].groupby("flask"):
            axo.plot(d.time_h, d.do_pct_as, color=CURATED_RAMP[2], lw=0.7)
        r = rep[rep.condition == c]
        for a_ in (ax, axp, axo):
            if data.VINESTOCK[c]["glucose_gl"] > 0:
                a_.axvspan(r.t_last_with_glucose_h.mean(), r.t_first_without_glucose_h.mean(), color=G200, lw=0, zorder=0)
                a_.axvline(r.t_glucose_depletion_h.mean(), color=G800, lw=0.7, ls=(0, (3, 2)))
            if data.VINESTOCK[c]["fructose_gl"] > 0:
                a_.axvline(r.t_depletion_h.mean(), color=G500, lw=0.7, ls=(0, (1, 1.5)))
        ax.set_ylim(-0.03, 0.85); axp.set_ylim(-0.035 * top, top); axo.set_ylim(0, 225)
        ax.set_xlim(0, np.ceil((on[on.condition == c].time_h.max() + 2) / 10) * 10)
        ax.tick_params(labelbottom=False); axp.tick_params(labelbottom=False)
        ax.set_title(FLASK[c], pad=3)
        axo.set_xlabel("Time (h)")
        if j == 0:
            ax.set_ylabel("Glycerol\n(g L$^{-1}$)")
            axp.set_ylabel("Arabitol\n(area µL$^{-1}$)")
            axo.set_ylabel("Dissolved O$_2$\n(% air sat.)")
    fig.legend([Line2D([], [], color=CMP["glycerol"], marker="o", ms=2.8), Line2D([], [], color=C16, marker="s", ms=2.6, mfc="white", mew=0.6),
                Line2D([], [], color=CURATED_RAMP[2], lw=0.7), Line2D([], [], color=G800, lw=0.7, ls=(0, (3, 2))),
                Patch(color=G200, lw=0), Line2D([], [], color=G500, lw=0.7, ls=(0, (1, 1.5)))],
               ["Glycerol", "Arabitol", "Dissolved O$_2$ (3 replicates)", "Glucose depleted", "Last sample with glucose to first without", "All sugars depleted"],
               loc="lower center", bbox_to_anchor=(0.52, 0.0), ncol=3, handlelength=1.8, columnspacing=1.6, labelspacing=0.3)
    save(fig, "fig4_products", ["flask_products", "flask_hplc"])


def figS_arabitol():
    """The 16.1 min peak against single standards: retention on the flask run and UV-to-RI ratio (results/peak_identity)."""
    tr = pd.read_csv(RES / "peak_identity" / "fig_traces.csv")
    uv = pd.read_csv(RES / "peak_identity" / "fig_uv.csv")
    order = ["Glc 15, 48 h", "Arabitol", "Xylitol", "Sorbitol", "Mannitol", "Malic acid"]
    shown = {"Glc 15, 48 h": "Glc 15, 48 h", "Arabitol": "Arabitol standard", "Xylitol": "Xylitol standard", "Sorbitol": "Sorbitol standard",
             "Mannitol": "Mannitol standard", "Malic acid": "Malic acid standard"}
    C16 = CMP["peak16"]
    fig = plt.figure(figsize=(139 * MM, 62 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.62, left=0.17, right=0.985, top=0.93, bottom=0.2)
    a, b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    panel(a, "a", x=-0.3, y=1.0); panel(b, "b", x=-0.62, y=1.0)
    n = len(order)
    for i, lab in enumerate(order):
        q = tr[tr.trace == lab]
        colour = C16 if lab.startswith("Glc") else (INK if lab == "Arabitol" else G500)
        a.plot(q.rt_on_flask_run, n - 1 - i + 0.8 * q.signal_scaled, color=colour, lw=0.9 if colour != G500 else 0.7)
    x0 = float(tr[tr.kind == "sample"].apex_on_flask_run.iloc[0])
    a.axvline(x0, color=G600, lw=0.6, ls=(0, (3, 2)), zorder=0)
    a.set_yticks([n - 1 - i + 0.1 for i in range(n)])
    a.set_yticklabels([shown[o] for o in order])
    for tick, o in zip(a.get_yticklabels(), order):
        tick.set_color(INK if o in ("Glc 15, 48 h", "Arabitol") else G700)
    a.tick_params(axis="y", length=0)
    a.spines["left"].set_visible(False)
    a.set_xlim(14.0, 17.4); a.set_ylim(-0.15, n)
    a.set_xlabel("Retention time (min)")
    ref = uv[uv.kind != "sample"].sort_values("uv210_ri", ascending=False)
    labels = list(ref.label) + ["16.1 min peak"]
    pos = np.arange(len(labels))[::-1]
    for p_, (_, r) in zip(pos[:-1], ref.iterrows()):
        b.plot([r.uv210_ri], [p_], "o", ms=3.4, color=INK if r.label == "Arabitol" else G500, mec="white", mew=0.5, zorder=3)
    s_ = uv[uv.kind == "sample"].uv210_ri.to_numpy()
    b.plot(s_, pos[-1] + np.linspace(-0.25, 0.25, len(s_)), "s", ms=2.6, color=C16, mfc="white", mew=0.6, ls="none", zorder=3)
    b.set_yticks(pos); b.set_yticklabels(labels)
    b.tick_params(axis="y", length=0)
    b.spines["left"].set_visible(False)
    b.set_xscale("symlog", linthresh=0.1, linscale=0.6)
    b.set_xlim(-0.25, 25)
    b.set_xticks([-0.1, 0, 0.1, 1, 10]); b.set_xticklabels(["\u22120.1", "0", "0.1", "1", "10"])
    b.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    b.grid(axis="x", color=G200, lw=0.5); b.set_axisbelow(True)
    b.set_xlabel("UV 210 nm / RI area")
    save(fig, "figS3_arabitol_identity", ["peak_identity"])


# ── Supplementary: yield sensitivity ──────────────────────────────────────────

def figS_sensitivity():
    df = pd.read_csv(RES / "yield_sensitivity" / "sensitivity.csv")
    df = df[df.sugar == "glucose"]
    rep = pd.read_csv(RES / "flask_products" / "replicates.csv")
    meas = {"Glc 15": rep[rep.condition == "Glc15"].yield_gg.mean(), "Fru 15": rep[rep.condition == "Fru15"].yield_gg.mean()}
    order = df[(df.model == "curated") & (df.change == "+25%")].assign(a=lambda d: d.yield_change_pct.abs()).sort_values("a").parameter.tolist()
    fig, axs = plt.subplots(1, 2, figsize=(FULL, 62 * MM), sharey=True, gridspec_kw={"wspace": 0.08, "left": 0.2, "right": 0.985, "top": 0.93, "bottom": 0.2})
    for ax, m, letter in zip(axs, MODELS, "ab"):
        panel(ax, letter, x=-0.04 if letter == "b" else -0.42)
        q = df[df.model == m]
        base = q.yield_base_gg.iloc[0]
        for i, par in enumerate(order):
            r = q[q.parameter == par].set_index("change").yield_change_pct
            ax.barh(i, r["+25%"], color=MODEL_C[m], height=0.6)
            ax.barh(i, r["-25%"], color=MODEL_C[m], height=0.6, alpha=0.45)
        for lab, y in meas.items():
            v = 100 * (y / base - 1)
            ax.axvline(v, color=MEAS, lw=0.7, ls=(0, (2, 1.5)))
            ax.text(v, len(order) - 0.35, lab, fontsize=5.6, ha="center", va="bottom")
        ax.axvline(0, color=G800, lw=0.6)
        ax.set_yticks(range(len(order)), [p_.replace("rowth-associated ATP", "rowth-associated maintenance") for p_ in order], fontsize=6.2)
        ax.set_title(MODEL_SHORT[m], pad=9)
        ax.set_xlim(-18, 6)
        ax.set_xlabel("Change in aerobic biomass yield on glucose (%)")
    fig.legend([Patch(color=G600), Patch(color=G600, alpha=0.45), Line2D([], [], color=MEAS, lw=0.7, ls=(0, (2, 1.5)))],
               ["parameter +25%", "parameter −25%", "change needed to reach the measured yield"], loc="lower center",
               bbox_to_anchor=(0.5, -0.04), ncol=3, handlelength=1.6, columnspacing=1.4)
    save(fig, "figS2_yield_sensitivity", ["yield_sensitivity", "flask_products"])


# ── Figure 5: pH without fitting ──────────────────────────────────────────────

def fig4():
    d = jload("ph_validation/series.json")
    S, A = d["summary_stage_a"], d["stage_a"]
    fig = plt.figure(figsize=(FULL, 108 * MM))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.05], hspace=0.6, wspace=0.42, left=0.065, right=0.985, top=0.96, bottom=0.1)
    ax = fig.add_subplot(gs[0, :2]); panel(ax, "a", x=-0.12)
    steps = [("v0", "Buffer\n+ acids"), ("v1", "+ NH$_4^+$\nuptake"),
             ("v2", "+ activity\ncorrection"), ("v3", "+ P, S\nuptake")]
    vals = [S[k]["rmse"] for k, _ in steps]
    ax.plot(range(4), vals, color=G400, lw=0.8, zorder=1)
    for i, (k, lab) in enumerate(steps):
        ax.scatter(i, vals[i], s=26, color=PHV[k], zorder=3, edgecolor="white", lw=0.5)
        ax.text(i, vals[i] * 1.32, f"{vals[i]:.3f}", ha="center", fontsize=6.3)
    ax.set_yscale("log"); ax.set_ylim(0.03, 2.0); ax.set_xlim(-0.4, 3.4)
    ax.set_xticks(range(4), [l for _, l in steps], fontsize=6)
    ax.set_ylabel("Mean pH RMSE (12 replicates)")
    ax = fig.add_subplot(gs[0, 2:]); panel(ax, "b", x=-0.12)
    for cond, fl in A.items():
        for f in fl.values():
            x, ph, v3 = (np.array(f[k], float) for k in ("x", "ph", "v3"))
            fr = np.clip((x - x[0]) / max(x[-1] - x[0], 1e-6), 0, 1)
            ax.scatter(fr, v3 - ph, s=2.2, color=MODEL_TINT["curated"], alpha=0.8, lw=0)
    allr = []
    for cond, fl in A.items():
        for f in fl.values():
            x, ph, v3 = (np.array(f[k], float) for k in ("x", "ph", "v3"))
            fr = np.clip((x - x[0]) / max(x[-1] - x[0], 1e-6), 0, 1); allr.append(np.c_[fr, v3 - ph])
    allr = np.vstack(allr); bins = np.linspace(0, 1, 11)
    mids = [(bins[i] + bins[i + 1]) / 2 for i in range(10)]
    mean = [np.nanmean(allr[(allr[:, 0] >= bins[i]) & (allr[:, 0] < bins[i + 1] + (i == 9)), 1]) for i in range(10)]
    ax.plot(mids, mean, color=MEAS, lw=1.0)
    ax.axhline(0, color=G600, lw=0.6)
    ax.set_ylim(-0.25, 0.25); ax.set_xlim(0, 1)
    ax.set_xlabel("Growth progress (fraction of final biomass)"); ax.set_ylabel("Predicted − measured pH\n(full model)")
    ax.legend([Line2D([], [], color=MEAS)], ["10% bin mean"], loc="lower right")
    for j, cond in enumerate(data.PRIMARY_VINESTOCK):
        ax = fig.add_subplot(gs[1, j]); panel(ax, "cdef"[j], x=-0.3)
        fl = A[cond]
        for f in fl.values():
            ax.plot(f["t"], f["ph"], color=G600, lw=0.6)
        num = lambda a: np.array([np.nan if y is None else y for y in a], float)
        tg = num(list(fl.values())[0]["t"])
        for v, ls in (("v0", (0, (2.5, 1.5))), ("v3", "-")):
            # mean of the three flask predictions on the first flask's time grid
            ys = []
            for fv in fl.values():
                t, y = num(fv["t"]), num(fv[v]); ok = ~np.isnan(t) & ~np.isnan(y)
                ys.append(np.interp(tg, t[ok], y[ok], left=np.nan, right=np.nan))
            ax.plot(tg, np.nanmean(ys, axis=0), color=PHV[v], lw=1.2 if v == "v3" else 0.9, ls=ls)
        ax.set_title(FLASK[cond], pad=3)
        ax.set_ylim(5.75, 6.8); ax.set_xlabel("Time (h)")
        if j == 0:
            ax.set_ylabel("Medium pH")
    fig.legend([Line2D([], [], color=G600, lw=0.6), Line2D([], [], color=PHV["v3"], lw=1.2),
                Line2D([], [], color=PHV["v0"], ls=(0, (2.5, 1.5)))],
               ["measured, online (3 replicates)", "predicted, full model (mean of replicates)", "predicted, buffer + acids only"],
               loc="lower center", bbox_to_anchor=(0.5, -0.05), ncol=3, handlelength=1.8, columnspacing=1.6)
    save(fig, "fig5_ph", ["ph_validation"])


# ── Tables ────────────────────────────────────────────────────────────────────

def tables(qc, st, gpr):
    """CSV versions of Table 1 (inputs), Table 2 (models) and Supplementary Table S1 (curation rules). The other
    supplementary tables are read directly from their result files (see the README)."""
    TAB.mkdir(parents=True, exist_ok=True)
    NM = {k: v[0] for k, v in names.MODEL.items()}          # plain-text names for CSV
    M = MODELS
    t1 = [["Reactions"] + [qc[m]["reactions"] for m in M], ["Metabolites"] + [qc[m]["metabolites"] for m in M],
          ["Genes"] + [qc[m]["genes"] for m in M],
          ["Exchange reactions open for secretion"] + [qc[m]["open_exports"] for m in M],
          ["Biomass weight (g mmol⁻¹)"] + [round(qc[m]["biomass_weight"], 2) for m in M],
          ["Gene associations whose gene EC matches the reaction EC"] + [gpr[m]["EC matches"] for m in M],
          ["Gene associations whose gene EC shares only the first three levels"] + [gpr[m]["related EC (same x.y.z)"] for m in M],
          ["Gene associations whose gene EC differs"] + [gpr[m]["different EC"] for m in M],
          ["Maximum ATP yield (mol ATP per mol glucose)"] + [round(qc[m]["max_atp_per_glucose"], 1) for m in M],
          ["Energy-generating cycle tests failed"] + [len(qc[m]["energy_leaks"]) for m in M]] + \
         [[f"Aerobic biomass yield on {s} (g g⁻¹)"] + [round(qc[m]["yield_gg"][s], 3) for m in M] for s in ("glucose", "fructose")]
    t1 = [[r[0]] + [f"{v:,}" if isinstance(v, int) else (f"{v:g}" if isinstance(v, float) else v) for v in r[1:]] for r in t1]
    pd.DataFrame(t1, columns=["Property"] + [NM[m] for m in M]).to_csv(TAB / "table2_models.csv", index=False)

    I = simulate.INPUTS
    meta = jload("ph_validation/series.json")["meta"]
    t2 = [
        ["Template model", "iJB1325 (A. niger), CBS 513.88 version", "literature", "Brandl et al. 2018, Additional file 3"],
        ["Biomass composition", "Fractions of growing mycelium (Supplementary Table S2), 1 g mmol⁻¹", "measured (protein, RNA), derived (others)",
         "Towersey et al. 1977; Trinci 1992; FDA GRN 000091; Pedersen et al. 1999"],
        ["Non-growth-associated maintenance", "1.9 mmol ATP gDW⁻¹ h⁻¹", "literature (A. niger)", "Andersen et al. 2008"],
        ["Growth-associated maintenance", "61 mmol ATP gDW⁻¹", "literature (A. niger)", "Andersen et al. 2008"],
        ["P/O ratio", "2.46 for NADH, 1.64 for succinate", "literature (template)", "12 and 8 protons pumped, 4.88 per cytosolic ATP"],
        ["Acid export cost", "1 ATP per negative charge at pH 7", "literature", "van Maris et al. 2004a, 2004b; Abbott et al. 2009"],
        ["Growth rate that sets the maximum sugar uptake", f"{I['mu_max_per_h'][0]:.2f} h⁻¹", I["mu_max_per_h"][1],
         "growth rate of the Quorn process, 0.17-0.20 h⁻¹, Wiebe 2002"],
        ["Sugar Km", f"{I['km_sugar_mM'][0]} mM", I["km_sugar_mM"][1],
         "between the glucose (0.025 mM) and fructose (4.5 mM) values of A. niger MstA, van Kuyk et al. 2004"],
        ["Maximum oxygen transfer rate", f"{data.VINESTOCK_OTR_MAX_MMOL_L_H:.1f} mmol L⁻¹ h⁻¹", "measured (kLa)",
         f"kLa {data.VINESTOCK_KLA_PER_H:g} h⁻¹, 200% air (Vinestock & Guo 2026)"],
        ["Initial biomass", "0.01 g L⁻¹ in the base settings; fitted per condition", "assumed; fitted",
         "time to half final biomass, which sets the lag time"],
        ["Sequential sugar uptake", f"other sugars blocked above {I['diauxie_threshold_mM'][0]} mM glucose", "assumed",
         "carbon catabolite repression, Adnan et al. 2018"],
        ["Initial sugars", "15 g L⁻¹ in total, as prepared", "measured (recipe)", "Vinestock & Guo 2026"],
        ["Medium", f"KH₂PO₄ {meta['medium_mM']['phosphate']:.0f} mM, NH₄Cl {meta['medium_mM']['NH4Cl']:.0f} mM", "measured (recipe)",
         "Vinestock & Guo 2026 (KH₂PO₄), Whittaker 2022 (NH₄Cl)"],
        ["Ionic strength", f"{meta['ionic_strength_M']:.3f} M, phosphate pKa₂′ {meta['pka2_apparent']:.2f}", "measured (recipe)",
         "calculated with the Davies equation and the charge balance"],
    ]
    pd.DataFrame(t2, columns=["Input", "Value", "Origin", "Source"]).to_csv(TAB / "table1_inputs.csv", index=False)

    tS1 = [["1", "Orthology filter: reactions kept only with orthology support (tier 1 RBH, tier 2 orthogroup, complexes kept when one subunit maps)", "Orthology", f"{st['draft_removed']} template reactions removed"],
           ["2", "Reaction recovery by EC number and gap-filling for biomass synthesis", "UniProt EC of the F. venenatum proteome; tasks of Thiele & Palsson 2010", f"+{st['ec_rescue']} by EC number, +{st['gapfill']} by gap-filling"],
           ["3", "Secretion restricted to a gene-supported transporter or literature support", "List of metabolites with citations", f"{st['exports_closed']} exchange reactions closed"],
           ["4", "Redox rules against cost-free NADH-NADPH and cytosol-mitochondria routes, and NAD(P) reaction directionality", "eQuilibrator |ln Γ| ≥ 7 (Noor et al. 2012)", f"{st['redox']} redox actions, {st['thermo_v1']} directionality changes"],
           ["5", "Biomass normalised to the composition of growing mycelium (Supplementary Table S2), NGAM 1.9", "Towersey et al. 1977; FDA GRN 000091; Andersen et al. 2008", "1 g mmol⁻¹"],
           ["6", "Thermodynamic directionality of all metabolic reactions", "eQuilibrator |ln Γ| ≥ 7 (Noor et al. 2012)", f"{st['thermo_v2']} directionality changes (incl. phosphoketolase)"],
           ["6", "Organic-acid export pays 1 ATP per negative charge", "van Maris et al. 2004a, 2004b; Abbott et al. 2009", f"{st['acid']} transports made import-only, {st['acid']} ATP-coupled export reactions added"],
           ["6", "Mitochondrial GTP/GDP carrier imports GTP", "Vozza et al. 2004", "1 reaction"]]
    pd.DataFrame(tS1, columns=["Stage", "Rule", "Basis", "Effect"]).to_csv(TAB / "tableS1_curation_rules.csv", index=False)


def only(names_: list[str]) -> None:
    """Redraw the named figure functions and keep the manifest of the others."""
    OUT.mkdir(parents=True, exist_ok=True)
    for name in names_:
        globals()[name]()
    path = OUT / "manifest.json"
    old = json.loads(path.read_text()) if path.exists() else {}
    path.write_text(json.dumps({**old, **MANIFEST}, indent=1))
    print("figures:", sorted(MANIFEST))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    qc = jload("model_versions_qc/qc.json")
    for m, f in (("draft", DRAFT), ("curated", CURATED)):
        qc[m].update(model_counts(f))
    st, _, gpr = fig1(qc)
    fig2(qc)
    fig3()
    fig_products()
    fig4()
    figS_sensitivity()
    figS_arabitol()
    tables(qc, st, gpr)
    (OUT / "manifest.json").write_text(json.dumps(MANIFEST, indent=1))
    print("figures:", sorted(MANIFEST))


if __name__ == "__main__":
    if len(sys.argv) > 1:
        only(sys.argv[1:])
    else:
        main()
