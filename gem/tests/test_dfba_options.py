import cobra
import pytest

from dfba import run_dfba_multi


def _toy():
    m = cobra.Model("toy")
    S_e, S, O2_e, O2, B = (cobra.Metabolite(i, compartment=c) for i, c in
                           [("S_e", "e"), ("S", "c"), ("O2_e", "e"), ("O2", "c"), ("BIOMASS", "c")])
    def add(rid, st, b=(0, 1000)):
        r = cobra.Reaction(rid); r.bounds = b; m.add_reactions([r]); r.add_metabolites(st); return r
    add("BOUNDARY_Se", {S_e: -1}, (-10, 1000))
    add("BOUNDARY_O2e", {O2_e: -1}, (-1000, 1000))
    add("T", {S_e: -1, O2_e: -1, S: 1, O2: 1})
    add("GROW", {S: -1, O2: -1, B: 1})
    add("DRAIN_Biomass", {B: -1})
    m.objective = "DRAIN_Biomass"
    return m


def test_o2_cap_fn_is_applied_each_step():
    m = _toy()
    seen = []
    def cap(t, x):
        seen.append(x)
        return 0.5
    sub = [{"name": "s", "rxn": m.reactions.BOUNDARY_Se, "mw": 1.0, "vmax": 10, "km": 0.1, "conc_mm": 100.0}]
    res = run_dfba_multi(m, m.reactions.DRAIN_Biomass, sub, 0.1, 1.0, 0.3, 0.1, [], o2_cap_fn=cap)
    assert len(seen) == len(res) and all(r["mu_fba"] <= 0.5 + 1e-9 for r in res)


def test_continue_on_infeasible_runs_to_end():
    m = _toy()
    m.reactions.GROW.lower_bound = 0.1          # forced flux -> infeasible once S is gone
    sub = [{"name": "s", "rxn": m.reactions.BOUNDARY_Se, "mw": 1.0, "vmax": 10, "km": 0.1, "conc_mm": 0.0}]
    res = run_dfba_multi(m, m.reactions.DRAIN_Biomass, sub, 0.1, 1.0, 0.5, 0.1, [],
                         continue_on_infeasible=True)
    assert len(res) >= 5 and all(r.get("infeasible") for r in res)


def _yield_model():
    """1 mmol substrate -> 1 g biomass: biomass gained must equal substrate used."""
    m = cobra.Model("y")
    S_e, S, B = (cobra.Metabolite(i, compartment=c) for i, c in [("S_e", "e"), ("S", "c"), ("BIOMASS", "c")])
    for rid, st, b in [("BOUNDARY_Se", {S_e: -1}, (-1000, 1000)), ("T", {S_e: -1, S: 1}, (0, 1000)),
                       ("GROW", {S: -1, B: 1}, (0, 1000)), ("DRAIN_Biomass", {B: -1}, (0, 1000))]:
        r = cobra.Reaction(rid); r.bounds = b; m.add_reactions([r]); r.add_metabolites(st)
    m.objective = "DRAIN_Biomass"
    return m


@pytest.mark.parametrize("dt", [0.5, 2.0])
def test_dfba_conserves_mass_for_any_step(dt):
    m = _yield_model()
    s0 = 10.0  # mM
    sub = [{"name": "s", "rxn": m.reactions.BOUNDARY_Se, "mw": 1.0, "vmax": 1.0, "km": 0.01, "conc_mm": s0}]
    res = run_dfba_multi(m, m.reactions.DRAIN_Biomass, sub, 0.1, 1.0, 40.0, dt, [], continue_on_infeasible=True)
    end = res[-1]
    gained = end["biomass"] - 0.1          # g
    used = s0 - end["c_BOUNDARY_Se"]        # mmol (1 L)
    assert used <= s0 + 1e-9
    assert gained == pytest.approx(used, rel=1e-6)
