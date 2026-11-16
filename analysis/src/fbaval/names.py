"""Display names used in every figure, table and caption.

Terms follow the reconstruction protocol (Thiele & Palsson 2010: draft reconstruction, curated model). Internal keys
(the ``model`` column of the result tables, result folders) stay short.
"""

# models: internal key -> (full name, short name)
MODEL = {
    "draft": ("Draft", "Draft"),
    "curated": ("Curated", "Curated"),
}
STAGE6 = "stage 6: thermodynamic directionality for all reactions and the energy cost of organic-acid export"

# cumulative pH chemistry: internal key -> (full name, short name)
PH = {
    "v0": ("Buffer + acids", "Buffer + acids"),
    "v1": ("+ NH₄⁺ uptake", "+ NH₄⁺ uptake"),
    "v2": ("+ activity correction", "+ activity corr."),
    "v3": ("+ P, S uptake (full model)", "+ P, S uptake"),
}

FLASK_COND = {"Glc15": "Glc 15", "Glc10Fru5": "Glc 10 + Fru 5", "Glc5Fru10": "Glc 5 + Fru 10", "Fru15": "Fru 15"}
CONC_NOTE = "Sugar concentrations in g L⁻¹ (Glc glucose, Fru fructose)."


def flask(cond: str) -> str:
    """'Glc 15' for a batch fermentation condition."""
    return FLASK_COND[cond]
