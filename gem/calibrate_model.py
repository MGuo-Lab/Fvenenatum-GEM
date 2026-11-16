"""Add a non-growth-associated maintenance (NGAM) ATP demand to the model.

iJB1325 / iMA871 only models growth-associated ATP costs (61 mmol ATP per gDW
biomass via the GROWTH reaction's GAM). Real cells also burn ATP for survival
processes that don't scale with growth (membrane potential, protein turnover,
hyphal extension). Without NGAM, FBA over-predicts biomass yield by
funnelling almost all carbon into biomass.

Default NGAM = 1.9 mmol ATP/gDW/h — Andersen et al. (2008), the original
iMA871 reconstruction paper, which estimated this from A. niger
glucose-limited chemostat data. F. venenatum-specific NGAM has not been
measured directly; the A. niger value is the closest published reference for
a filamentous Ascomycete with similar physiology, and it is cited here as a
literature default rather than as an organism-specific parameter.

GAM (61 mmol ATP/gDW) is already encoded in the iJB1325 GROWTH reaction
stoichiometry — this script only adds the non-growth-associated component.

The reaction added is:

    NGAM:  ATP + H2O --> ADP + H + PI    (lb = ub = ngam_flux)

with bounds equal so the flux is fixed regardless of the FBA objective.

References:
- Andersen MR, Nielsen ML, Nielsen J. *Metabolic model integration of the
  bibliome, genome, metabolome and reactome of Aspergillus niger.*
  Mol Syst Biol. 2008;4:178. (iMA871 paper; mATP = 1.9 mmol ATP/gDW/h)
- Brandl J, Andersen MR. *A community-driven reconstruction of the
  Aspergillus niger metabolic network.* Fungal Biol Biotechnol. 2018;5:16.
  (iJB1325 paper; the template used in this pipeline)

Usage::

    python calibrate_model.py --input build/model.xml --ngam 1.9 \
        --output build/model_calibrated.xml
"""

from __future__ import annotations

import argparse
import os
import sys

import cobra


# ── NGAM addition ────────────────────────────────────────────────────────────

def add_ngam(model: cobra.Model, flux: float, rxn_id: str = "NGAM") -> None:
    """Add (or update) a forced ATP hydrolysis NGAM reaction."""

    if rxn_id in model.reactions:
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = flux
        rxn.upper_bound = flux
        print(f"Updated existing {rxn_id}: lb=ub={flux}")
        return

    # Resolve the cytosolic ATP pool. iMA871 calls it "ATP" (cytoplasm) and
    # "ATPm" (mitochondria); cytosolic NGAM is the standard choice.
    needed = ["ATP", "H2O", "ADP", "PI", "H"]
    missing = [m for m in needed if m not in model.metabolites]
    if missing:
        sys.exit(f"ERROR: NGAM substrates not found in model: {missing}")

    rxn = cobra.Reaction(rxn_id)
    rxn.name = "Non-growth-associated maintenance ATP hydrolysis"
    rxn.lower_bound = flux
    rxn.upper_bound = flux
    rxn.add_metabolites({
        model.metabolites.ATP:  -1.0,
        model.metabolites.H2O:  -1.0,
        model.metabolites.ADP:  +1.0,
        model.metabolites.PI:   +1.0,
        model.metabolites.H:    +1.0,
    })
    model.add_reactions([rxn])
    print(f"Added {rxn_id}: ATP + H2O --> ADP + H + PI    (lb=ub={flux})")


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--input",  default="build/model.xml")
    parser.add_argument("--output", default="build/model_calibrated.xml")
    parser.add_argument("--ngam",   type=float, default=1.9,
                        help="NGAM ATP hydrolysis flux (mmol ATP/gDW/h). Default 1.9 from "
                             "Andersen et al. 2008 (A. niger iMA871, glucose-limited chemostat). "
                             "F. venenatum-specific value has not been measured.")
    parser.add_argument("--probe-medium", action="store_true",
                        help="Print growth & yield with the standard glucose minimal "
                             "medium so the calibration impact is visible.")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f"ERROR: input model not found: {args.input}")

    print(f"Loading: {args.input}")
    m = cobra.io.read_sbml_model(args.input)
    print(f"  {len(m.reactions)} reactions, {len(m.metabolites)} metabolites, "
          f"{len(m.genes)} genes")

    add_ngam(m, args.ngam)

    if args.probe_medium:
        # Quick FBA on glucose minimal to show effect of calibration
        print("\nProbe: glucose minimal medium (lb=-10)...")
        for r in m.exchanges:
            r.lower_bound = 0.0
        for rid, lb in [
            ("BOUNDARY_GLCe", -10), ("BOUNDARY_NH4e", -20),
            ("BOUNDARY_PIe",  -20), ("BOUNDARY_SLFe",  -20),
            ("BOUNDARY_O2e",  -20), ("BOUNDARY_H2Oe",  -1000),
        ]:
            if rid in m.reactions:
                m.reactions.get_by_id(rid).lower_bound = lb
        sol = m.optimize()
        if sol.status == "optimal":
            mu = sol.objective_value
            glc = -sol.fluxes.get("BOUNDARY_GLCe", 0.0)
            print(f"  μ = {mu:.4f} 1/h")
            print(f"  glucose uptake = {glc:.4f} mmol/gDW/h")
            if glc > 1e-6:
                # Approximate biomass-on-substrate yield in g/g
                # (rough: μ in g/(gDW·h), glc in mmol/(gDW·h), MW=180.156 g/mol)
                y_gg = mu / (glc * 180.156 / 1000.0)
                print(f"  Y_biomass/glucose ≈ {y_gg:.3f} g/g")

    # Restore the model to its bounds before save (probe-medium modified bounds in-place)
    if args.probe_medium:
        print(f"\nReloading {args.input} for save (preserving original bounds)...")
        m = cobra.io.read_sbml_model(args.input)
        add_ngam(m, args.ngam)

    cobra.io.write_sbml_model(m, args.output)
    print(f"\nCalibrated model saved to: {args.output}")


if __name__ == "__main__":
    main()
