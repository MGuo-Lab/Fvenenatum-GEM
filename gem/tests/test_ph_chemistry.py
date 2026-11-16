import math

import cobra
import pytest

from dfba_kinetics import (CardinalPH, apparent_pkas, calibrate_net_cation, make_acid,
                           step_pH_state)
from dfba import run_dfba_multi


def _hh_phosphate(pH0, total_mM, acid_mM, pka2=7.20):
    """Henderson-Hasselbalch around pKa2 after adding `acid_mM` of strong acid."""
    r0 = 10 ** (pH0 - pka2)
    b0 = total_mM * r0 / (1 + r0)
    a0 = total_mM - b0
    return pka2 + math.log10((b0 - acid_mM) / (a0 + acid_mM))


def test_no_assimilation_keeps_initial_pH():
    cation = calibrate_net_cation(6.6, [make_acid("phosphate", 200.0)])
    assert step_pH_state({}, 0.0, 200.0, cation) == pytest.approx(6.6, abs=1e-3)


def test_ammonium_assimilation_acidifies_like_strong_acid():
    # NH4+ is a strong cation below pH 8 (pKa 9.25): each NH4+ built into biomass
    # removes one positive charge from the medium.
    cation = calibrate_net_cation(6.6, [make_acid("phosphate", 200.0)])
    pH = step_pH_state({}, 0.0, 200.0, cation, ammonium_assimilated_mM=30.0)
    assert pH == pytest.approx(_hh_phosphate(6.6, 200.0, 30.0), abs=0.02)


def test_davies_lowers_phosphate_pka2_at_high_ionic_strength():
    assert apparent_pkas((2.15, 7.20, 12.35), 0.0) == pytest.approx((2.15, 7.20, 12.35))
    pk = apparent_pkas((2.15, 7.20, 12.35), ionic_strength_M=0.39, temp_c=28.0)
    assert pk[1] == pytest.approx(6.79, abs=0.03)


def test_ionic_strength_is_threaded_through_the_solver():
    # pKa2' (about 6.8) is nearer the start pH than 7.2, so the buffer is stronger and
    # the same acid load moves the pH less.
    c0 = calibrate_net_cation(6.6, [make_acid("phosphate", 200.0)])
    c1 = calibrate_net_cation(6.6, [make_acid("phosphate", 200.0, ionic_strength_M=0.39)])
    d0 = 6.6 - step_pH_state({}, 0.0, 200.0, c0, ammonium_assimilated_mM=30.0)
    d1 = 6.6 - step_pH_state({}, 0.0, 200.0, c1, ammonium_assimilated_mM=30.0,
                             ionic_strength_M=0.39)
    assert 0 < d1 < d0


def _model(acetate_per_g=0.0):
    """1 mmol S + 1 mmol NH4 -> 1 g biomass (+ optional acetate)."""
    m = cobra.Model("n")
    S_e, S, N_e, N, A, A_e, B = (cobra.Metabolite(i, compartment=c) for i, c in
                                 [("S_e", "e"), ("S", "c"), ("NH4_e", "e"), ("NH4", "c"),
                                  ("AC", "c"), ("AC_e", "e"), ("BIOMASS", "c")])
    grow = {S: -1, N: -1, B: 1}
    if acetate_per_g:
        grow[A] = acetate_per_g
    for rid, st, b in [("BOUNDARY_Se", {S_e: -1}, (-1000, 1000)),
                       ("BOUNDARY_NH4e", {N_e: -1}, (-1000, 1000)),
                       ("BOUNDARY_ACe", {A_e: -1}, (0, 1000)),
                       ("TS", {S_e: -1, S: 1}, (0, 1000)), ("TN", {N_e: -1, N: 1}, (0, 1000)),
                       ("TA", {A: -1, A_e: 1}, (0, 1000)),
                       ("GROW", grow, (0, 1000)), ("DRAIN_Biomass", {B: -1}, (0, 1000))]:
        r = cobra.Reaction(rid); r.bounds = b; m.add_reactions([r]); r.add_metabolites(st)
    m.objective = "DRAIN_Biomass"
    return m


def test_dfba_tracks_ammonium_for_pH():
    m = _model()
    sub = [{"name": "s", "rxn": m.reactions.BOUNDARY_Se, "mw": 1.0, "vmax": 1.0, "km": 0.01,
            "conc_mm": 20.0}]
    res = run_dfba_multi(m, m.reactions.DRAIN_Biomass, sub, 0.1, 1.0, 10.0, 0.1, [],
                         pH_params=CardinalPH(pH_min=2.0, pH_opt=6.6, pH_max=9.0),
                         initial_pH=6.6, medium_phosphate_mM=200.0,
                         continue_on_infeasible=True)
    last = res[-1]
    assert last["nh4_assimilated_mM"] > 1.0
    assert last["pH"] == pytest.approx(_hh_phosphate(6.6, 200.0, last["nh4_assimilated_mM"]),
                                       abs=0.03)


def test_dfba_counts_acid_byproducts_without_overflow():
    m = _model(acetate_per_g=2.0)
    m.reactions.BOUNDARY_NH4e.bounds = (-1000, 1000)
    sub = [{"name": "s", "rxn": m.reactions.BOUNDARY_Se, "mw": 1.0, "vmax": 1.0, "km": 0.01,
            "conc_mm": 20.0}]
    kw = dict(pH_params=CardinalPH(pH_min=2.0, pH_opt=6.6, pH_max=9.0), initial_pH=6.6,
              medium_phosphate_mM=200.0, nh4_rxn_id=None, continue_on_infeasible=True)
    without = run_dfba_multi(m, m.reactions.DRAIN_Biomass, sub, 0.1, 1.0, 10.0, 0.1,
                             ["BOUNDARY_ACe"], **kw)
    sub[0]["conc_mm"] = 20.0
    with_acid = run_dfba_multi(m, m.reactions.DRAIN_Biomass, sub, 0.1, 1.0, 10.0, 0.1,
                               ["BOUNDARY_ACe"], acid_keys={"BOUNDARY_ACe": "acetate"}, **kw)
    assert without[-1]["pH"] == pytest.approx(6.6, abs=1e-3)
    assert with_acid[-1]["pH"] < 6.5


def test_infeasible_steps_still_report_medium_pH():
    m = _model()
    sub = [{"name": "s", "rxn": m.reactions.BOUNDARY_Se, "mw": 1.0, "vmax": 1.0, "km": 0.01,
            "conc_mm": 2.0}]
    m.reactions.GROW.lower_bound = 0.01                 # infeasible once S runs out
    res = run_dfba_multi(m, m.reactions.DRAIN_Biomass, sub, 0.1, 1.0, 30.0, 0.5, [],
                         pH_params=CardinalPH(pH_min=2.0, pH_opt=6.6, pH_max=9.0),
                         initial_pH=6.6, medium_phosphate_mM=200.0,
                         continue_on_infeasible=True)
    dead = [r for r in res if r.get("infeasible")]
    assert dead and all(r["pH"] == r["pH"] and r["pH"] < 6.6 for r in dead)


def _spill_model():
    """S + NH4 -> 1 g biomass, or S -> W (an exportable sink); NH4 acidifies."""
    m = _model()
    S = m.metabolites.get_by_id("S")
    W, W_e = cobra.Metabolite("W", compartment="c"), cobra.Metabolite("W_e", compartment="e")
    for rid, st, b in [("SPILL", {S: -1, W: 1}, (0, 1000)), ("TW", {W: -1, W_e: 1}, (0, 1000)),
                       ("BOUNDARY_We", {W_e: -1}, (0, 1000))]:
        r = cobra.Reaction(rid); r.bounds = b; m.add_reactions([r]); r.add_metabolites(st)
    return m


def _run_spill(coupling):
    m = _spill_model()
    sub = [{"name": "s", "rxn": m.reactions.BOUNDARY_Se, "mw": 1.0, "vmax": 1.0, "km": 0.01,
            "conc_mm": 40.0}]
    return run_dfba_multi(m, m.reactions.DRAIN_Biomass, sub, 0.5, 1.0, 12.0, 0.1, ["BOUNDARY_We"],
                          pH_params=CardinalPH(pH_min=5.0, pH_opt=6.6, pH_max=9.0),
                          initial_pH=6.6, medium_phosphate_mM=20.0, pH_coupling=coupling,
                          continue_on_infeasible=True)


def test_growth_cap_coupling_keeps_uptake_and_spills_carbon():
    # Acidification slows growth in both modes, but only "growth_cap" keeps taking
    # up substrate at the kinetic rate and sends the surplus to the export.
    scale, cap = _run_spill("scale_all"), _run_spill("growth_cap")
    s_left = lambda res: res[-1]["c_BOUNDARY_Se"]
    assert s_left(cap) < s_left(scale)
    assert cap[-1]["c_BOUNDARY_We"] > 0 and scale[-1]["c_BOUNDARY_We"] == pytest.approx(0, abs=1e-6)
    assert min(r["gamma_pH"] for r in cap) < 0.9           # the pH brake was active


def test_default_coupling_is_unchanged():
    m = _model()
    sub = [{"name": "s", "rxn": m.reactions.BOUNDARY_Se, "mw": 1.0, "vmax": 1.0, "km": 0.01,
            "conc_mm": 20.0}]
    kw = dict(pH_params=CardinalPH(pH_min=2.0, pH_opt=6.6, pH_max=9.0), initial_pH=6.6,
              medium_phosphate_mM=200.0, continue_on_infeasible=True)
    a = run_dfba_multi(m, m.reactions.DRAIN_Biomass, sub, 0.1, 1.0, 5.0, 0.1, [], **kw)
    sub[0]["conc_mm"] = 20.0
    b = run_dfba_multi(m, m.reactions.DRAIN_Biomass, sub, 0.1, 1.0, 5.0, 0.1, [],
                       pH_coupling="scale_all", **kw)
    assert a[-1]["biomass"] == pytest.approx(b[-1]["biomass"])


def test_growth_cap_keeps_sugar_uptake_after_nitrogen_runs_out():
    # NH4 is an on-demand substrate: once it is gone growth stops, but the sugar
    # is still taken up at its kinetic rate and leaves as the spill product.
    m = _spill_model()
    sub = [{"name": "s", "rxn": m.reactions.BOUNDARY_Se, "mw": 1.0, "vmax": 1.0, "km": 0.01,
            "conc_mm": 30.0},
           {"name": "nh4", "rxn": m.reactions.BOUNDARY_NH4e, "mw": 1.0, "vmax": 100.0,
            "km": 0.01, "conc_mm": 5.0, "on_demand": True}]
    res = run_dfba_multi(m, m.reactions.DRAIN_Biomass, sub, 1.0, 1.0, 20.0, 0.1,
                         ["BOUNDARY_We"], pH_coupling="growth_cap",
                         continue_on_infeasible=True)
    last = res[-1]
    assert last["biomass"] == pytest.approx(1.0 + 5.0, rel=1e-3)   # 1 g per mmol N
    assert last["c_BOUNDARY_Se"] < 0.1
    assert last["c_BOUNDARY_We"] == pytest.approx(30.0 - 5.0, rel=1e-2)
    assert not any(r.get("coupling_fallback") for r in res)


def test_glucose_first_rule_does_not_block_on_demand_nutrients():
    # The diauxic rule holds back secondary sugars, not the nitrogen source.
    m = _spill_model()
    sub = [{"name": "s", "rxn": m.reactions.BOUNDARY_Se, "mw": 1.0, "vmax": 1.0, "km": 0.01,
            "conc_mm": 10.0},
           {"name": "nh4", "rxn": m.reactions.BOUNDARY_NH4e, "mw": 1.0, "vmax": 100.0,
            "km": 0.01, "conc_mm": 50.0, "on_demand": True}]
    res = run_dfba_multi(m, m.reactions.DRAIN_Biomass, sub, 1.0, 1.0, 5.0, 0.1, [],
                         diauxic_primary_id="BOUNDARY_Se", pH_coupling="growth_cap",
                         continue_on_infeasible=True)
    assert res[-1]["biomass"] > 2.0


def _sink_model():
    """Spill model + sulfate needed for growth (0.01 per g) and an inorganic
    electron sink: S + SLF -> CO2 + sulfite, which avoids exporting organic W."""
    m = _spill_model()
    S, B = m.metabolites.S, m.metabolites.BIOMASS
    SL_e, SL, SO3, SO3_e, C, C_e = (cobra.Metabolite(i, compartment=c) for i, c in
                                    [("SLF_e", "e"), ("SLF", "c"), ("SO3", "c"), ("SO3_e", "e"),
                                     ("CO2", "c"), ("CO2_e", "e")])
    C.formula, C_e.formula = "CO2", "CO2"
    for rid, st, b in [("BOUNDARY_SLFe", {SL_e: -1}, (-1000, 1000)), ("TSL", {SL_e: -1, SL: 1}, (0, 1000)),
                       ("DUMP", {S: -1, SL: -1, C: 1, SO3: 1}, (0, 1000)),
                       ("TSO3", {SO3: -1, SO3_e: 1}, (0, 1000)), ("BOUNDARY_SO3e", {SO3_e: -1}, (0, 1000)),
                       ("TC", {C: -1, C_e: 1}, (0, 1000)), ("BOUNDARY_CO2e", {C_e: -1}, (0, 1000))]:
        r = cobra.Reaction(rid); r.bounds = b; m.add_reactions([r]); r.add_metabolites(st)
    m.reactions.GROW.add_metabolites({SL: -0.01})
    return m


def test_growth_cap_holds_mineral_nutrients_to_growth_needs():
    # Only carbon is uncoupled from growth: sulfate (like NH4 and phosphate) is taken up
    # in proportion to growth, so it cannot serve as an electron sink for surplus carbon.
    m = _sink_model()
    sub = [{"name": "s", "rxn": m.reactions.BOUNDARY_Se, "mw": 1.0, "vmax": 1.0, "km": 0.01,
            "conc_mm": 40.0}]
    res = run_dfba_multi(m, m.reactions.DRAIN_Biomass, sub, 0.5, 1.0, 12.0, 0.1, ["BOUNDARY_We"],
                         pH_params=CardinalPH(pH_min=5.0, pH_opt=6.6, pH_max=9.0),
                         initial_pH=6.6, medium_phosphate_mM=20.0, pH_coupling="growth_cap",
                         continue_on_infeasible=True)
    grown = res[-1]["biomass"] - 0.5
    assert res[-1]["c_BOUNDARY_We"] > 1.0                 # surplus leaves as the organic spill
    assert res[-1]["so4_assimilated_mM"] == pytest.approx(0.01 * grown, rel=1e-3)
