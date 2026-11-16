"""Locations inside this repository.

The package is installed in editable mode from ``analysis/`` (see the README), so the repository root is three
levels above this file. ``FBAVAL_RESULTS`` redirects the results folder, for example to compare a rerun with the
committed results. ``FBAVAL_RAW_DIR`` points to the raw instrument files, which are not part of the repository: the
HPLC data files (Shimadzu LabSolutions ``.lcd``) under ``lcd/`` and the PreSens exports under ``presens/``.
"""

from __future__ import annotations

import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
MODEL_DIR = REPO / "model"
CURATION_DIR = MODEL_DIR / "curation"
GEM_DIR = REPO / "gem"
GEM_INPUTS = GEM_DIR / "inputs"
DATA_DIR = REPO / "data"
FLASK_DIR = DATA_DIR / "flask"
RESULTS = Path(os.environ.get("FBAVAL_RESULTS", str(REPO / "results")))
FIGURES = REPO / "figures"
RAW_DIR = Path(os.environ.get("FBAVAL_RAW_DIR", str(REPO / "raw")))

DRAFT = MODEL_DIR / "Draft.xml"
CURATED = MODEL_DIR / "iJL1143.xml"
TEMPLATE = GEM_DIR / "raw_data" / "MOESM3_cbs51388.xml"    # downloaded by gem/download_data.py


def model_file(key: str) -> Path:
    """SBML file of a model by its key in the result tables: ``draft`` or ``curated`` (iJL1143)."""
    return {"draft": DRAFT, "curated": CURATED}[key]
