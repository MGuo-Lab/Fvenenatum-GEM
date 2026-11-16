"""MEMOTE snapshot scores for the Draft and Curated models.

Outputs: results/memote_ijl1143/{draft.json, iJL1143.json, scores.csv, provenance.json}

Usage (from analysis/):  python scripts/run_memote.py
"""

from __future__ import annotations

import json

import memote
import pandas as pd
from memote.suite import api

from fbaval import simulate
from fbaval.provenance import write_provenance

from fbaval.paths import CURATED, DRAFT, RESULTS

OUT = RESULTS / "memote_ijl1143"
MODELS = {"draft": DRAFT, "iJL1143": CURATED}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for key, fname in MODELS.items():
        path = str(fname)
        model = simulate.load(path)
        code, result = api.test_model(model, sbml_version=(3, 1), results=True, pytest_args=["--tb=no", "-q", "-p", "no:cacheprovider"], solver_timeout=120)
        report = json.loads(api.snapshot_report(result, html=False))
        (OUT / f"{key}.json").write_text(json.dumps(report))
        score = report["score"]
        rows.append({"model": key, "section": "total", "score": score["total_score"]})
        rows += [{"model": key, "section": s["section"], "score": s["score"]} for s in score["sections"]]
        print(key, "exit", code, "total", round(score["total_score"], 3))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "scores.csv", index=False)
    write_provenance(OUT, params={"tool": f"memote {memote.__version__}", "models": {k: v.name for k, v in MODELS.items()}, "solver_timeout_s": 120, "fitted_parameters": "none"})
    print(df.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
