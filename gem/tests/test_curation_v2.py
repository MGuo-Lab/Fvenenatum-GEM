"""Curation v2 rules: acid-export energetics and the thermodynamic direction
rule widened from NAD(P) reactions to all metabolic reactions."""

import cobra
import pytest

import curate_model as cm


def _rxn(model, rid, stoich, bounds=(0, 1000)):
    r = cobra.Reaction(rid)
    r.bounds = bounds
    model.add_reactions([r])
    r.add_metabolites(stoich)
    return r


@pytest.fixture
def acid_model():
    """S -> AC + 2 ATP; AC leaves through a reversible transport (uptake as written)."""
    m = cobra.Model("acid")
    mets = {k: cobra.Metabolite(k, compartment=c) for k, c in
            [("S_e", "e"), ("S", "c"), ("AC", "c"), ("ACe", "e"), ("ATP", "c"),
             ("ADP", "c"), ("PI", "c"), ("H2O", "c")]}
    g = mets.__getitem__
    _rxn(m, "BOUNDARY_Se", {g("S_e"): -1}, (-1, 1000))
    _rxn(m, "T_S", {g("S_e"): -1, g("S"): 1})
    _rxn(m, "FERM", {g("S"): -1, g("AC"): 1, g("ADP"): -2, g("PI"): -2, g("ATP"): 2, g("H2O"): 2})
    _rxn(m, "r1069", {g("ACe"): -1, g("AC"): 1}, (-1000, 1000))
    _rxn(m, "BOUNDARY_ACe", {g("ACe"): -1})
    _rxn(m, "ATPM", {g("ATP"): -1, g("H2O"): -1, g("ADP"): 1, g("PI"): 1})
    m.objective = "ATPM"
    return m


def test_acid_export_pays_one_atp_per_charge(acid_model):
    assert acid_model.slim_optimize() == pytest.approx(2.0)
    log = cm.CurationLog()
    changed = cm.charge_acid_export(acid_model, {"AC": 1}, log)
    assert changed == ["r1069"]
    assert acid_model.slim_optimize() == pytest.approx(1.0)     # 2 ATP made, 1 spent exporting
    t = acid_model.reactions.r1069
    assert t.bounds == (0.0, 1000.0)                            # uptake only, still free
    exp = acid_model.reactions.get_by_id("r1069_exp")
    assert exp.metabolites[acid_model.metabolites.ATP] == -1
    assert any(r["reaction_id"] == "r1069" and r["stage"] == "6_v2_acid_export" for r in log.rows)


def test_acid_export_handles_transport_written_outward(acid_model):
    t = acid_model.reactions.r1069
    t.add_metabolites({acid_model.metabolites.ACe: 2, acid_model.metabolites.AC: -2})  # AC -> ACe
    assert cm.charge_acid_export(acid_model, {"AC": 2}, cm.CurationLog()) == ["r1069"]
    assert t.bounds == (-1000.0, 0.0)                           # only the inward direction left
    assert acid_model.reactions.r1069_exp.metabolites[acid_model.metabolites.ATP] == -2


class _Q:
    def __init__(self, v, e):
        self.value, self.error = type("M", (), {"magnitude": v})(), type("M", (), {"magnitude": e})()


class _FakeCC:
    """eQuilibrator stand-in: every reaction strongly favours the forward direction."""
    def parse_reaction_formula(self, s):
        return type("Rx", (), {"is_balanced": lambda self, ignore_atoms=(): True})()

    def ln_reversibility_index(self, rx):
        return _Q(-12.0, 0.3)


def _annotated_model():
    m = cobra.Model("thermo")
    A, B, N, NH = (cobra.Metabolite(i, compartment="c") for i in ("A", "B", "NAD", "NADH"))
    for met, kegg in ((A, "C00231"), (B, "C00227"), (N, "C00003"), (NH, "C00004")):
        met.annotation["kegg.compound"] = kegg
    Bx = cobra.Metabolite("BIOMASS", compartment="c")
    _rxn(m, "EX_A", {A: -1}, (-10, 1000))
    _rxn(m, "PK", {A: -1, B: 1}, (-1000, 1000))                  # no NAD(P)
    _rxn(m, "DH", {A: -1, N: -1, B: 1, NH: 1}, (-1000, 1000))     # NAD-linked
    _rxn(m, "RECYC", {NH: -1, N: 1})
    _rxn(m, "GROW", {B: -1, Bx: 1})
    _rxn(m, "DRAIN_Biomass", {Bx: -1})
    m.objective = "DRAIN_Biomass"
    return m


def test_thermo_rule_scope():
    m = _annotated_model()
    med = {"EX_A": -10}
    cm.thermodynamic_directions(m, med, cm.CurationLog(), cc=_FakeCC())
    assert m.reactions.PK.bounds == (-1000, 1000)               # v1: NAD(P) reactions only
    assert m.reactions.DH.bounds == (0, 1000)
    m = _annotated_model()
    cm.thermodynamic_directions(m, med, cm.CurationLog(), cc=_FakeCC(), redox_only=False)
    assert m.reactions.PK.bounds == (0, 1000)                    # v2: every metabolic reaction
