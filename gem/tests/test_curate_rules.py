import cobra
import pytest

import curate_model as cm


def _met(mid, comp, formula=None):
    m = cobra.Metabolite(mid, compartment=comp, formula=formula)
    return m


def _rxn(model, rid, stoich, bounds=(0, 1000), gpr=""):
    r = cobra.Reaction(rid)
    r.bounds = bounds
    model.add_reactions([r])
    r.add_metabolites(stoich)
    r.gene_reaction_rule = gpr
    return r


@pytest.fixture
def toy():
    """S_e -> S -> BIOMASS (+ W byproduct that must leave); X, Y exportable."""
    m = cobra.Model("toy")
    S_e, S, W, W_e, X, X_e, Y, Y_e, B = (
        _met("S_e", "e"), _met("S", "c"), _met("W", "c"), _met("W_e", "e"),
        _met("X", "c"), _met("X_e", "e"), _met("Y", "c"), _met("Y_e", "e"),
        _met("BIOMASS", "c"))
    _rxn(m, "BOUNDARY_Se", {S_e: -1}, (-10, 1000))
    _rxn(m, "T_S", {S_e: -1, S: 1}, gpr="g1")
    _rxn(m, "R_grow", {S: -1, B: 1, W: 1}, gpr="g2")
    _rxn(m, "DRAIN_Biomass", {B: -1})
    _rxn(m, "T_W", {W: -1, W_e: 1})                 # gene-less transporter
    _rxn(m, "BOUNDARY_We", {W_e: -1})
    _rxn(m, "T_X", {X: -1, X_e: 1}, gpr="g3")       # gene-backed transporter
    _rxn(m, "BOUNDARY_Xe", {X_e: -1})
    _rxn(m, "T_Y", {Y: -1, Y_e: 1})                 # gene-less, whitelisted
    _rxn(m, "BOUNDARY_Ye", {Y_e: -1})
    m.objective = "DRAIN_Biomass"
    return m


MED = {"BOUNDARY_Se": -10}


def test_secretion_rule_closes_only_unsupported_and_reopens_required(toy):
    log = cm.CurationLog()
    closed = cm.apply_secretion_rule(toy, {"BOUNDARY_Ye": "diffusion; ref"}, MED, log)
    assert toy.reactions.BOUNDARY_Xe.upper_bound > 0          # gene-backed
    assert toy.reactions.BOUNDARY_Ye.upper_bound > 0          # whitelisted
    assert toy.reactions.BOUNDARY_We.upper_bound > 0          # reopened: needed for growth
    assert closed == []
    assert any(r["action"] == "export_reopened" and r["reaction_id"] == "BOUNDARY_We"
               for r in log.rows)


def test_secretion_rule_closes_unneeded_geneless_export(toy):
    W_e = toy.metabolites.W_e
    Z, Z_e = _met("Z", "c"), _met("Z_e", "e")
    _rxn(toy, "T_Z", {Z: -1, Z_e: 1})
    _rxn(toy, "BOUNDARY_Ze", {Z_e: -1})
    closed = cm.apply_secretion_rule(toy, {}, MED, cm.CurationLog())
    assert "BOUNDARY_Ze" in closed
    assert toy.reactions.BOUNDARY_Ze.upper_bound == 0
    assert W_e is toy.metabolites.W_e


def test_needs_based_gapfill_adds_only_required(toy):
    template = toy.copy()
    toy.remove_reactions([toy.reactions.T_S])
    extra = _rxn(template, "R_extra", {template.metabolites.S: -1, template.metabolites.X: 1})
    log = cm.CurationLog()
    def biomass(m):
        with m:
            cm.set_medium(m, MED)
            return cm.growth(m) >= cm.GROWTH_THRESHOLD
    added = cm.needs_based_gapfill(toy, template, {"biomass": biomass}, log)
    assert added == ["T_S"]
    assert "R_extra" not in toy.reactions and extra.id == "R_extra"
    assert log.rows[0]["confidence"] == 1


def test_ec_rescue_uses_proteome_ec():
    template = cobra.Model("t")
    a, b = _met("A", "c"), _met("B", "c")
    r1 = _rxn(template, "R1", {a: -1, b: 1}, gpr="tg1")
    r1.annotation["ec-code"] = "1.1.1.1"
    r2 = _rxn(template, "R2", {b: -1, a: 1}, gpr="tg2")
    r2.annotation["ec-code"] = "9.9.9.9"
    model = cobra.Model("m")
    added = cm.ec_rescue(model, template, {"1.1.1.1": ["P1", "P2"]}, cm.CurationLog())
    assert added == ["R1"]
    assert model.reactions.R1.gene_reaction_rule == "P1 or P2"


def test_block_pyridine_nucleotide_transport():
    m = cobra.Model("m")
    nadh, nadhm = _met("NADH", "c"), _met("NADHm", "m")
    _rxn(m, "r_nadh_t", {nadh: -1, nadhm: 1})
    _rxn(m, "r_nadh_t_gene", {nadh: -1, nadhm: 1}, gpr="g1")
    blocked = cm.block_pyridine_nucleotide_transport(m, cm.CurationLog())
    assert blocked == ["r_nadh_t"]
    assert m.reactions.r_nadh_t.bounds == (0, 0)


def test_mitochondrial_adh_made_reductive():
    m = cobra.Model("m")
    eth, nad, acal, nadh = (_met("ETHm", "m"), _met("NADm", "m"),
                            _met("ACALm", "m"), _met("NADHm", "m"))
    r = _rxn(m, "r118a", {eth: -1, nad: -1, acal: 1, nadh: 1}, (-1000, 1000))
    r.annotation["ec-code"] = "1.1.1.1"
    cm.mitochondrial_adh_reductive(m, cm.CurationLog())
    assert r.bounds == (-1000, 0)   # only acetaldehyde + NADHm -> ethanol + NADm


def test_biomass_weight_expands_pseudo_metabolites():
    m = cobra.Model("m")
    glc = _met("GLC", "c", "C6H12O6")
    h2o = _met("H2O", "c", "H2O")
    poly = _met("POLY", "c")
    bio = _met("BIOMASS", "c")
    _rxn(m, "polym", {glc: -1, h2o: 1, poly: 1})
    _rxn(m, "GROWTH", {poly: -2, bio: 1})
    w = cm.biomass_weight(m)
    assert w == pytest.approx(2 * (180.156 - 18.015) / 1000, rel=1e-3)


def _pyridine_model():
    m = cobra.Model("m")
    mets = {i: _met(i, "c") for i in ("GL", "GLYN", "NAD", "NADH", "NADP", "NADPH")}
    return m, mets


def test_remove_gene_less_cofactor_variant():
    m, x = _pyridine_model()
    _rxn(m, "r122", {x["GL"]: -1, x["NAD"]: -1, x["GLYN"]: 1, x["NADH"]: 1}, (-1000, 1000))
    _rxn(m, "r127", {x["GL"]: -1, x["NADP"]: -1, x["GLYN"]: 1, x["NADPH"]: 1}, (-1000, 1000),
         gpr="F1")
    removed = cm.remove_cofactor_duplicates(m, cm.CurationLog())
    assert removed == ["r122"]
    assert "r127" in m.reactions


def test_keep_cofactor_variants_when_both_gene_backed():
    m, x = _pyridine_model()
    _rxn(m, "a", {x["GL"]: -1, x["NAD"]: -1, x["GLYN"]: 1, x["NADH"]: 1}, gpr="F1")
    _rxn(m, "b", {x["GL"]: -1, x["NADP"]: -1, x["GLYN"]: 1, x["NADPH"]: 1}, gpr="F2")
    assert cm.remove_cofactor_duplicates(m, cm.CurationLog()) == []


def test_transhydrogenase_needs_ec_annotated_gene():
    m, x = _pyridine_model()
    _rxn(m, "th_bad", {x["NADPH"]: -1, x["NAD"]: -1, x["NADP"]: 1, x["NADH"]: 1},
         (-1000, 1000), gpr="MFS1")
    _rxn(m, "th_ok", {x["NADPH"]: -1, x["NAD"]: -1, x["NADP"]: 1, x["NADH"]: 1},
         (-1000, 1000), gpr="TH1")
    removed = cm.remove_unsupported_transhydrogenases(m, {"1.6.1.2": ["TH1"]}, cm.CurationLog())
    assert removed == ["th_bad"]


def test_orient_cofactor_pairs():
    m, x = _pyridine_model()
    nad = _rxn(m, "m_nad", {x["GL"]: -1, x["NAD"]: -1, x["GLYN"]: 1, x["NADH"]: 1},
               (-1000, 1000), gpr="F1")
    nadp = _rxn(m, "m_nadp", {x["GL"]: -1, x["NADP"]: -1, x["GLYN"]: 1, x["NADPH"]: 1},
                (-1000, 1000), gpr="F2")
    changed = cm.orient_cofactor_pairs(m, {}, cm.CurationLog())
    assert sorted(changed) == ["m_nad", "m_nadp"]
    assert nad.bounds == (0, 1000)      # GL + NAD -> GLYN + NADH only
    assert nadp.bounds == (-1000, 0)    # GLYN + NADPH -> GL + NADP only


def test_normalize_biomass_sets_mass_fractions():
    def build(c_a, c_b):
        m = cobra.Model("m")
        glc, h2o = _met("GLC", "c", "C6H12O6"), _met("H2O", "c", "H2O")
        a, b, bio = _met("A", "c"), _met("B", "c"), _met("BIOMASS", "c")
        _rxn(m, "make_A", {glc: -1, h2o: 1, a: 1})          # A = 0.162 g/mmol
        _rxn(m, "make_B", {glc: -2, h2o: 2, b: 1})          # B = 0.324 g/mmol
        _rxn(m, "GROWTH", {a: -c_a, b: -c_b, bio: 1})
        return m
    model, template = build(1.0, 1.0), build(1.0, 1.0)
    new = cm.normalize_biomass(model, template, {"A": 30, "B": 70}, cm.CurationLog())
    assert new["A"] * 0.16214 == pytest.approx(0.3, rel=1e-3)
    assert cm.biomass_weight(model) == pytest.approx(1.0, rel=1e-3)


def test_orient_cofactor_pairs_skips_irreversible_pairs():
    m, x = _pyridine_model()
    a = _rxn(m, "ndh", {x["GL"]: -1, x["NAD"]: -1, x["GLYN"]: 1, x["NADH"]: 1}, (0, 1000), gpr="F1")
    b = _rxn(m, "ndph", {x["GL"]: -1, x["NADP"]: -1, x["GLYN"]: 1, x["NADPH"]: 1}, (0, 1000), gpr="F2")
    assert cm.orient_cofactor_pairs(m, {}, cm.CurationLog()) == []
    assert a.bounds == (0, 1000) and b.bounds == (0, 1000)


def test_formula_weight_approximates_r_groups_as_c18_1():
    pe = _met("PE", "c", "C7H12NO8PR2")
    assert cm.formula_weight(pe) == pytest.approx(744.0, abs=1.5)   # dioleoyl-PE
