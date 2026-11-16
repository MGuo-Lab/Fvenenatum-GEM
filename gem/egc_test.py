"""Energy-generating-cycle (EGC) / energy-leak test for the F. venenatum GEM.

Standard Fritzemeier et al. 2017 (PLOS Comp Biol) / Thiele & Palsson 2010
diagnostic: close ALL nutrient uptakes, then for each energy currency add a
"dissipation" reaction that discharges it (ATP -> ADP+Pi, NAD(P)H -> NAD(P),
QH2 -> Q, proton gradient -> cytosol) and maximise it. With nothing coming in,
a thermodynamically correct network can sustain ZERO flux through any of these.
Any positive flux = the model charges an energy metabolite from nothing = a leak
that inflates growth/ATP yield (field average ~25% biomass inflation).

Usage:  python egc_test.py [model.xml]
Read-only: operates on a temporary copy of each reaction via `with model:`.
Import ``find_energy_leaks`` to run the same tests from other scripts.
"""
from __future__ import annotations

import sys
import cobra
from cobra import Reaction

MODEL = sys.argv[1] if len(sys.argv) > 1 else "../model/iJL1143.xml"  # CLI only
TOL = 1e-6  # mmol/gDW/h below which we call it "no flux"


def by_name(model, name, comp="c"):
    """Return the metabolite whose name == `name` (case-insensitive) in compartment `comp`."""
    name = name.strip().lower()
    cands = [m for m in model.metabolites if (m.name or "").strip().lower() == name]
    incomp = [m for m in cands if m.compartment == comp]
    return (incomp or cands or [None])[0]


def add_dissipation(model, mid, stoich):
    """Add a temporary dissipation reaction (objective) and return (max_flux, solution)."""
    rxn = Reaction("EGC_" + mid)
    rxn.lower_bound = 0.0
    rxn.upper_bound = 1000.0
    rxn.add_metabolites(stoich)
    model.add_reactions([rxn])
    model.objective = rxn
    sol = model.optimize()
    flux = sol.objective_value if sol.status == "optimal" else float("nan")
    return flux, sol


def find_energy_leaks(model: cobra.Model, extra_tests=None) -> list[tuple[str, float]]:
    """Return [(label, max_flux)] for every energy-currency dissipation test.

    Works on a context-managed copy of ``model`` (no lasting changes): closes
    all boundary uptakes, relaxes forced fluxes (e.g. NGAM) to 0, then
    maximises each dissipation reaction. ``extra_tests`` is an optional list
    of (label, {met_id: coeff}) dissipations to test as well. NaN = skipped.
    """
    results: list[tuple[str, float]] = []
    with model:
        for r in model.reactions:
            if r.boundary and r.lower_bound < 0:
                r.lower_bound = 0.0
            if r.lower_bound > 0:
                r.lower_bound = 0.0

        ids = {x.id for x in model.metabolites}

        def m(mid):
            return model.metabolites.get_by_id(mid) if mid in ids else None

        h2o = by_name(model, "H2O") or by_name(model, "Water")
        pi = m("PI")
        h_c = m("H")

        for tp, dp in [("ATP", "ADP"), ("GTP", "GDP"), ("CTP", "CDP"),
                       ("UTP", "UDP"), ("ITP", "IDP")]:
            a, b = m(tp), m(dp)
            if a and b and h2o and pi and h_c:
                with model:
                    flux, _ = add_dissipation(model, tp, {a: -1, h2o: -1, b: 1, pi: 1, h_c: 1})
                results.append((f"{tp} -> {dp} + Pi", flux))

        for red, ox in [("NADH", "NAD"), ("NADPH", "NADP"), ("FADH2p", "FADp"),
                        ("QH2m", "Qm"), ("NADHm", "NADm"), ("NADPHm", "NADPm"),
                        ("FMNH2", "FMN")]:
            a, b = m(red), m(ox)
            if not (a and b):
                continue
            nH = a.elements.get("H", 0) - b.elements.get("H", 0)
            other_ok = all(a.elements.get(e, 0) == b.elements.get(e, 0)
                           for e in set(a.elements) | set(b.elements) if e != "H")
            proton = m("Hm") if a.compartment == "m" else h_c
            if not other_ok or nH < 0 or proton is None:
                results.append((f"{red} -> {ox} (SKIP: non-clean couple)", float("nan")))
                continue
            with model:
                flux, _ = add_dissipation(model, red, {a: -1, b: 1, proton: nH})
            results.append((f"{red} -> {ox} + {nH}H", flux))

        # Dissipate the gradient downhill: protons are pumped out of the matrix
        # (H_POm -> H_PO, e.g. Complex I) and flow back through ATP synthase
        # (r245: H_PO -> H_POm), so the high-energy side is H_PO / H.
        for src, dst in [("H", "Hm"), ("H_PO", "H_POm")]:
            a, b = m(src), m(dst)
            if a and b:
                with model:
                    flux, _ = add_dissipation(model, src, {a: -1, b: 1})
                results.append((f"{src} -> {dst} (PMF)", flux))

        for label, stoich in (extra_tests or []):
            mets = {m(k): v for k, v in stoich.items()}
            if None in mets:
                results.append((f"{label} (SKIP: metabolite missing)", float("nan")))
                continue
            with model:
                flux, _ = add_dissipation(model, "X_" + label.replace(" ", "_")[:40], mets)
            results.append((label, flux))
    return results


def main():
    print(f"loading {MODEL}")
    model = cobra.io.read_sbml_model(MODEL)
    print(f"  {len(model.reactions)} reactions, {len(model.metabolites)} metabolites\n")
    results = find_energy_leaks(model)

    print("=" * 64)
    print(f"{'dissipation reaction':<34}{'max flux':>14}   verdict")
    print("-" * 64)
    n_leak = 0
    for label, flux in results:
        if flux != flux:  # nan
            verdict = "skipped"
        elif abs(flux) > TOL:
            verdict = "*** LEAK ***"
            n_leak += 1
        else:
            verdict = "clean (0)"
        fs = "  n/a" if flux != flux else f"{flux:>14.4g}"
        print(f"{label:<34}{fs}   {verdict}")
    print("=" * 64)
    print(f"\nRESULT: {n_leak} energy currency(ies) leak. "
          f"{'MODEL HAS EGCs.' if n_leak else 'No energy-generating cycles detected.'}")


if __name__ == "__main__":
    main()
