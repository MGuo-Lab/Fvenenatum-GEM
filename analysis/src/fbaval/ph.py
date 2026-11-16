"""pH of the Vinestock RHM shake flasks, predicted from medium chemistry.

Nothing here is fitted. The medium composition comes from the recipe
(Whittaker 2022 RHM, as modified in Vinestock & Guo Table 2). The amount
of N, P and S built into each gram of biomass comes from the biomass
reaction of iJL1143. The strong-ion balance is calibrated once, at
t = 0, to the measured initial pH. The chemistry itself is
`dfba_kinetics` in gem/ (charge balance, pKa table, Davies correction).

Model versions add one physical effect at a time, in an order fixed
before the comparison was run:

    v0  phosphate buffer + measured organic acids
    v1  + NH4+ assimilation (one strong cation removed per N in biomass)
    v2  + ionic-strength (Davies) apparent pKas at 28 °C
    v3  + phosphate and sulfate assimilation
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import dfba_kinetics as dk

# Flask medium, mM in the 50 mL broth: 46 mL RHM made with 30 g/L KH2PO4.
DILUTION = 46.0 / 50.0
PHOSPHATE_MM = 30.0 * DILUTION / 136.086 * 1000.0      # KH2PO4
NH4CL_MM = 4.4 * DILUTION / 53.491 * 1000.0
K2SO4_MM = 0.3 * DILUTION / 174.26 * 1000.0
MGSO4_MM = 0.25 * DILUTION / 246.47 * 1000.0
TEMP_C = 28.0


@dataclass(frozen=True)
class Medium:
    """Salts in the broth (mM) and the culture temperature."""
    phosphate: float        # KH2PO4
    nh4cl: float
    k2so4: float
    mgso4: float
    temp_c: float = 28.0


VINESTOCK_MEDIUM = Medium(PHOSPHATE_MM, NH4CL_MM, K2SO4_MM, MGSO4_MM, TEMP_C)


@dataclass(frozen=True)
class Version:
    ammonium: bool
    davies: bool
    p_and_s: bool


VERSIONS = {
    "v0": Version(ammonium=False, davies=False, p_and_s=False),
    "v1": Version(ammonium=True, davies=False, p_and_s=False),
    "v2": Version(ammonium=True, davies=True, p_and_s=False),
    "v3": Version(ammonium=True, davies=True, p_and_s=True),
}


@dataclass(frozen=True)
class Composition:
    """mmol of N, P and S assimilated per g biomass."""
    n: float
    p: float
    s: float


def biomass_composition(model) -> Composition:
    """N, P and S per gram from the model's own uptake at its aerobic optimum
    (fixed by the biomass reaction, so independent of growth rate)."""
    import curate_model as cm
    with model:
        medium = dict(cm.MINIMAL_GLUCOSE)
        cm.set_medium(model, medium)
        sol = model.optimize()
        mu = sol.objective_value
        f = sol.fluxes
    return Composition(-f["BOUNDARY_NH4e"] / mu, -f["BOUNDARY_PIe"] / mu, -f["BOUNDARY_SLFe"] / mu)


def ionic_strength(pH0: float, iterations: int = 20,
                   medium: Medium = VINESTOCK_MEDIUM) -> float:
    """Ionic strength (M) of the fresh medium at pH0. Na+ (from the
    NaOH used to set the pH) closes the charge balance; the phosphate
    split uses Davies pKas at the current estimate, iterated to a fixed point."""
    I = 0.3
    pi, nh4, mg = medium.phosphate, medium.nh4cl, medium.mgso4
    k = 2.0 * medium.k2so4
    so4 = medium.k2so4 + mg
    for _ in range(iterations):
        acid = dk.make_acid("phosphate", pi, I, medium.temp_c)
        z = -acid.avg_charge(pH0)                       # mean negative charge per phosphate
        h = 10.0 ** (-pH0)
        f = [1.0]
        kp = 1.0
        for pk in acid.pKas:
            kp *= 10.0 ** (-pk)
            f.append(kp / h ** len(f))
        frac = np.array(f) / sum(f)                     # H3PO4, H2PO4-, HPO4 2-, PO4 3-
        cations = pi + k + nh4 + 2 * mg                 # K+, NH4+, Mg2+ charge
        anions = pi * z + nh4 + 2 * so4                 # phosphate, Cl-, SO4 2-
        na = max(anions - cations, 0.0)
        sum_cz2 = (pi + k + nh4 + na + 4 * mg + nh4 + 4 * so4
                   + pi * sum(frac[i] * i * i for i in range(4)))
        I = 0.5 * sum_cz2 / 1000.0
    return I


def predict_ph(biomass_gl, acids_mM: dict[str, np.ndarray], initial_pH: float,
               version: Version, comp: Composition,
               medium: Medium = VINESTOCK_MEDIUM) -> np.ndarray:
    """pH at each sample. `biomass_gl` is biomass formed since t0 (g/L);
    `acids_mM` maps PKA_TABLE names to measured concentrations (mM)."""
    x = np.asarray(biomass_gl, dtype=float)
    I = ionic_strength(initial_pH, medium=medium) if version.davies else 0.0
    base = [dk.make_acid("phosphate", medium.phosphate, I, medium.temp_c)]
    cation = dk.calibrate_net_cation(initial_pH, base)
    out = np.empty_like(x)
    for i, xi in enumerate(x):
        acids = {k: float(v[i]) for k, v in acids_mM.items() if v[i] > 0}
        out[i] = dk.step_pH_state(
            acids, 0.0, medium.phosphate, cation,
            ammonium_assimilated_mM=comp.n * xi if version.ammonium else 0.0,
            phosphate_assimilated_mM=comp.p * xi if version.p_and_s else 0.0,
            sulfate_assimilated_mM=comp.s * xi if version.p_and_s else 0.0,
            ionic_strength_M=I, temp_c=medium.temp_c)
    return out
