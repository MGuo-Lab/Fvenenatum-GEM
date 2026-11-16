"""Draft model, the comparator of iJL1143: the template transferred to F. venenatum without curation.

The CBS 513.88 version of iJB1325 with its GPRs remapped to F. venenatum
through the ortholog table, keeping every reaction (unmapped ones keep a
stripped GPR), plus only the two settings the curated model also uses:
the biomass composition of growing mycelium at 1 g/mmol (BIOMASS_COMPOSITION in
curate_model.py) and NGAM 1.9 mmol ATP/gDW/h.
A comparison with the curated model therefore isolates the curation rules.

Usage (from gem/):
    python build_draft_baseline.py            # -> model/Draft.xml
"""

from __future__ import annotations

import argparse
import io
import os
from contextlib import redirect_stdout

import cobra

import build_draft_model as bdm
import curate_model as cm
from calibrate_model import add_ngam


HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    p.add_argument("--template", default=os.path.join(HERE, "raw_data", "MOESM3_cbs51388.xml"))
    p.add_argument("--orthologs", default=os.path.join(HERE, "inputs", "orthologs_cbs513.tsv"))
    p.add_argument("--ngam", type=float, default=1.9)
    p.add_argument("--output", default=os.path.join(HERE, "..", "model", "Draft.xml"))
    args = p.parse_args()

    template = cobra.io.read_sbml_model(args.template)
    with redirect_stdout(io.StringIO()):
        draft = bdm.build_draft_model(template, bdm.load_ortholog_mapping(args.orthologs),
                                      keep_unmapped=True)
        cm.normalize_biomass(draft, template, cm.BIOMASS_COMPOSITION, cm.CurationLog())
        add_ngam(draft, args.ngam)
    cm.set_medium(draft, cm.MINIMAL_GLUCOSE)
    draft.id = "fven_a3_5_draft"
    draft.name = "Fusarium venenatum A3/5 draft: iJB1325 (CBS 513.88) transferred, no curation"
    print(cm.summarize(draft, "draft baseline", cm.MINIMAL_GLUCOSE))
    print(f"  biomass weight {cm.biomass_weight(draft):.3f} g/mmol; max ATP/glucose {cm.max_atp_yield(draft):.1f}")
    cobra.io.write_sbml_model(draft, args.output)
    print("Saved", args.output)


if __name__ == "__main__":
    main()
