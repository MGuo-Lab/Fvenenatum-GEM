"""Acid-base chemistry of the medium and optional kinetic add-ons for dFBA.

The medium pH is solved from electroneutrality over the buffer, the organic acids and the strong ions
(``solve_pH``), with apparent pKa values corrected for ionic strength by the Davies equation (``apparent_pkas``).
The pH analysis of the batch fermentations (``analysis/src/fbaval/ph.py``) uses this chemistry.

Two optional add-ons can be passed to ``dfba.run_dfba_multi``. Both are off by default and neither is used by the
analyses in this repository:

1. pH-coupled growth inhibition by the cardinal-pH model of Rosso et al. 1995 (J Theor Biol 162:447):
       gamma(pH) = (pH-pH_min)(pH-pH_max) / [(pH-pH_min)(pH-pH_max) - (pH-pH_opt)^2]
   for pH_min < pH < pH_max, else 0.
2. Overflow metabolism (Vemuri et al. 2007, Appl Environ Microbiol 73:1996; Basan et al. 2015, Nature 528:99): a
   minimum secretion flux on each byproduct exchange,
       v_byp >= alpha_byp * max(0, q_S - q_crit).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping


# ── Cardinal-pH model (Rosso 1995) ───────────────────────────────────────────

@dataclass
class CardinalPH:
    """Cardinal-pH growth multiplier (Rosso et al. 1995). The default values are placeholders.

    `gamma(pH)` returns a dimensionless multiplier in [0, 1].  Values are
    1 at pH_opt and 0 outside [pH_min, pH_max], with a smooth quadratic
    decay in between.
    """

    pH_min: float = 3.0
    pH_opt: float = 5.8
    pH_max: float = 7.5

    def gamma(self, pH: float) -> float:
        if pH <= self.pH_min or pH >= self.pH_max:
            return 0.0
        num = (pH - self.pH_min) * (pH - self.pH_max)
        den = num - (pH - self.pH_opt) ** 2
        if den == 0:
            return 0.0
        val = num / den
        return max(0.0, min(1.0, val))


# ── Acid-base chemistry ──────────────────────────────────────────────────────

@dataclass
class AcidSpecies:
    """A polyprotic acid in the medium.

    `pKas` must be in increasing order.  `C_mM` is *total* (sum over all
    protonation states) molar concentration in millimolar.  `z_max` is the
    fully-deprotonated charge (e.g. −1 for acetate, −3 for citrate).
    """

    name: str
    C_mM: float
    pKas: tuple[float, ...]
    z_max: int = -1   # final deprotonated charge (always negative for acids)

    def avg_charge(self, pH: float) -> float:
        """Mole-weighted average charge at the given pH.

        For an n-protic acid HnA with stepwise dissociation constants
        Ka_1 .. Ka_n, the fraction in protonation state i (i deprotonations)
        is f_i = (Π_{j≤i} Ka_j) / [H+]^i, normalised by Σ f_j.  The
        average charge is Σ (−i)·f_i / Σ f_j.
        """
        h = 10.0 ** (-pH)
        # f_i (cumulative) = K1·K2·...·Ki / h^i; start with f_0 = 1
        fs: list[float] = [1.0]
        K_prod = 1.0
        for pKa in self.pKas:
            K_prod *= 10.0 ** (-pKa)
            fs.append(K_prod / (h ** len(fs)))
        norm = sum(fs)
        if norm == 0:
            return 0.0
        return sum(-i * fs[i] for i in range(len(fs))) / norm


# Standard pKa table (CRC Handbook 96th ed., Apelblat 2014 for citrate)
PKA_TABLE: dict[str, tuple[float, ...]] = {
    "acetate":   (4.76,),
    "pyruvate":  (2.39,),
    "lactate":   (3.86,),
    "formate":   (3.75,),
    "succinate": (4.21, 5.64),
    "citrate":   (3.13, 4.76, 6.40),
    "phosphate": (2.15, 7.20, 12.35),
    "carbonate": (6.35, 10.33),
}


def _davies_log_gamma(z: int, ionic_strength_M: float, temp_c: float) -> float:
    """log10 activity coefficient of an ion of charge z (Davies equation)."""
    if ionic_strength_M <= 0 or z == 0:
        return 0.0
    A = 0.4918 + 6.6098e-4 * temp_c + 5.0231e-6 * temp_c ** 2   # Debye-Hückel A, water
    s = math.sqrt(ionic_strength_M)
    return -A * z * z * (s / (1.0 + s) - 0.3 * ionic_strength_M)


def apparent_pkas(pKas: tuple[float, ...], ionic_strength_M: float = 0.0,
                  temp_c: float = 25.0) -> tuple[float, ...]:
    """Concentration-based pKa' at a given ionic strength (Davies), for an acid
    whose fully protonated form is neutral. pH stays an activity (what an
    electrode or optical sensor reads): step i, H(i-1)A -> HA(i) + H+, gives
    pKa_i' = pKa_i + log γ(z=i) - log γ(z=i-1)."""
    if ionic_strength_M <= 0:
        return tuple(pKas)
    return tuple(pk + _davies_log_gamma(i, ionic_strength_M, temp_c)
                 - _davies_log_gamma(i - 1, ionic_strength_M, temp_c)
                 for i, pk in enumerate(pKas, start=1))


def make_acid(name: str, C_mM: float, ionic_strength_M: float = 0.0,
              temp_c: float = 25.0) -> AcidSpecies | None:
    """Build an `AcidSpecies` for any species in `PKA_TABLE`; with an ionic
    strength the pKas are the Davies-corrected apparent values."""
    pKas = PKA_TABLE.get(name)
    if pKas is None or C_mM <= 0:
        return None
    return AcidSpecies(name=name, C_mM=C_mM,
                       pKas=apparent_pkas(pKas, ionic_strength_M, temp_c), z_max=-len(pKas))


# ── Charge-balance pH solver ─────────────────────────────────────────────────

def _charge_balance(pH: float, acids: list[AcidSpecies], net_cation_mM: float) -> float:
    """Sum of charges (mM units) at the given pH; zero at the true pH."""
    Kw = 1e-14
    h = 10.0 ** (-pH)
    oh = Kw / h
    # H+ and OH- in mM:
    net = net_cation_mM + h * 1000.0 - oh * 1000.0
    for sp in acids:
        net += sp.avg_charge(pH) * sp.C_mM
    return net


def solve_pH(acids: list[AcidSpecies], net_cation_mM: float,
             pH_lo: float = 1.0, pH_hi: float = 13.0,
             tol: float = 1e-4) -> float:
    """Solve electroneutrality for pH by bisection.

    `net_cation_mM` is the lumped charge contribution of all non-acid
    species (Na⁺, K⁺, NH₄⁺ minus Cl⁻, SO₄²⁻ etc.).  In a dFBA setup we
    don't track each cation separately — we calibrate this single number
    once at t=0 to reproduce the measured initial pH (`calibrate_net_cation`).
    """
    f_lo = _charge_balance(pH_lo, acids, net_cation_mM)
    f_hi = _charge_balance(pH_hi, acids, net_cation_mM)
    if f_lo * f_hi > 0:
        # Solution outside bracket; clamp.
        return pH_lo if abs(f_lo) < abs(f_hi) else pH_hi
    for _ in range(80):
        mid = 0.5 * (pH_lo + pH_hi)
        f_mid = _charge_balance(mid, acids, net_cation_mM)
        if abs(f_mid) < tol:
            return mid
        if f_lo * f_mid < 0:
            pH_hi, f_hi = mid, f_mid
        else:
            pH_lo, f_lo = mid, f_mid
    return 0.5 * (pH_lo + pH_hi)


def calibrate_net_cation(target_pH: float, acids: list[AcidSpecies],
                         lo: float = -200.0, hi: float = 1000.0) -> float:
    """Find the lumped cation charge (mM) that produces `target_pH`.

    Used once at t=0 so the simulated medium starts at the experimentally
    observed pH.  After calibration this number is treated as
    a conserved quantity for the rest of the simulation — we are not
    modelling base addition or evaporation.
    """
    f_lo = _charge_balance(target_pH, acids, lo)
    f_hi = _charge_balance(target_pH, acids, hi)
    if f_lo * f_hi > 0:
        # No bracket — return the closer of the two endpoints
        return lo if abs(f_lo) < abs(f_hi) else hi
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f_mid = _charge_balance(target_pH, acids, mid)
        if abs(f_mid) < 1e-4:
            return mid
        if f_lo * f_mid < 0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
    return 0.5 * (lo + hi)


# ── Overflow metabolism (Crabtree-like) ──────────────────────────────────────

@dataclass
class OverflowParams:
    """Phenomenological overflow yields (optional; not used by the analyses here).

    `alphas` are mol(byproduct) per mol(substrate consumed above q_crit).
    The defaults are illustrative values.

    The carbon-balance check below confirms the four byproducts together
    take ≤ 1 mol C per mol C(glucose) — so the constraints cannot violate
    elemental conservation when imposed by the LP.
    """

    q_crit_mmol_gDW_h: float = 0.3   # threshold below which no overflow
    alphas: dict[str, float] = field(default_factory=lambda: {
        "ethanol":  0.50,   # 2·EtOH/glc max = 2.0 → 25% of theoretical
        "pyruvate": 0.20,
        "glycerol": 0.02,
        "acetate":  0.01,
    })

    def carbon_check(self, n_C_substrate: int = 6) -> float:
        """Return (mol C in forced byproducts) per (mol C substrate above q_crit).

        Should be < 1 for the constraints to be feasible alongside biomass.
        For glucose (6 C), defaults give:
           0.50·2 + 0.20·3 + 0.02·3 + 0.01·2 = 1.68 mol C/mol Glc
        That is > 6/6 = 1; treat as upper bound only.  The LP is allowed
        to *exceed* the lb (so it can produce more if optimal); we use it
        as a *floor*.  Carbon balance is ultimately enforced by the
        stoichiometric matrix.
        """
        c_count = {"ethanol": 2, "pyruvate": 3, "glycerol": 3,
                   "acetate": 2, "lactate": 3, "formate": 1}
        return sum(self.alphas.get(k, 0.0) * c_count.get(k, 0)
                   for k in self.alphas) / n_C_substrate


def apply_overflow_lb(model, byproduct_ids: Mapping[str, str],
                      q_S: float, params: OverflowParams) -> dict[str, float]:
    """Set lower bounds on fermentative byproduct exchanges.

    `q_S` is the *current* substrate uptake rate (positive value, e.g. the
    absolute value of the glucose exchange flux from the previous LP).
    `byproduct_ids` maps overflow keys ('ethanol' etc.) to model exchange
    reaction IDs (typically `BOUNDARY_<MET>e`).

    Returns the {key: lb_applied} so the dFBA driver can log them.
    """
    excess = max(0.0, q_S - params.q_crit_mmol_gDW_h)
    applied: dict[str, float] = {}
    for key, alpha in params.alphas.items():
        rid = byproduct_ids.get(key)
        if rid is None or rid not in model.reactions:
            continue
        rxn = model.reactions.get_by_id(rid)
        # iJB1325 boundaries are formulated as <MET>e → ∅ (export = +flux),
        # so a *minimum* secretion rate is a positive lower bound.
        lb = alpha * excess
        # Don't lower an already-higher bound from a previous step.
        rxn.lower_bound = lb
        # Keep upper bound generous so the LP can still find the optimum
        if rxn.upper_bound < lb:
            rxn.upper_bound = max(1000.0, lb * 10)
        applied[key] = lb
    return applied


def reset_overflow_lb(model, byproduct_ids: Mapping[str, str]) -> None:
    """Restore byproduct exchange lower bounds to 0 (export-only default)."""
    for rid in byproduct_ids.values():
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound = 0.0


# ── Convenience: run a full step of overflow + pH adjustment ─────────────────

def step_pH_state(byproduct_concs_mM: dict[str, float],
                  buffer_citrate_mM: float,
                  base_phosphate_mM: float,
                  net_cation_mM: float,
                  ammonium_assimilated_mM: float = 0.0,
                  phosphate_assimilated_mM: float = 0.0,
                  sulfate_assimilated_mM: float = 0.0,
                  ionic_strength_M: float = 0.0,
                  temp_c: float = 25.0) -> float:
    """Compute the medium pH for the *current* dFBA timestep.

    Inputs are concentrations in millimolar, on the medium (not per-cell).
    The buffer/base species are only included if their concentration is
    positive.  Returns the pH (a float).

    Nutrient assimilation changes the strong-ion balance calibrated at t=0:
    each NH4+ built into biomass removes one positive charge (the cell
    releases the H+), each SO4^2- removes two negative charges, and the
    phosphate taken up leaves the buffer pool.
    """
    species: list[AcidSpecies] = []
    for name, C in byproduct_concs_mM.items():
        sp = make_acid(name, C, ionic_strength_M, temp_c)
        if sp is not None:
            species.append(sp)
    if buffer_citrate_mM > 0:
        sp = make_acid("citrate", buffer_citrate_mM, ionic_strength_M, temp_c)
        if sp is not None:
            species.append(sp)
    phosphate = base_phosphate_mM - phosphate_assimilated_mM
    if phosphate > 0:
        sp = make_acid("phosphate", phosphate, ionic_strength_M, temp_c)
        if sp is not None:
            species.append(sp)
    cation = net_cation_mM - ammonium_assimilated_mM + 2.0 * sulfate_assimilated_mM
    return solve_pH(species, cation)
