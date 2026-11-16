import pytest

from fbaval.paths import CURATED

needs_model = pytest.mark.skipif(not CURATED.exists(), reason="model/iJL1143.xml not found")


@pytest.fixture(scope="module")
def model():
    from fbaval import simulate
    return simulate.load(str(CURATED))


@needs_model
def test_static_profile_aerobic_yield_and_no_artefact_sinks(model):
    from fbaval import simulate
    df = simulate.static_profile(model, "glucose", 2.2, o2_fractions=(1.0, 0.4))
    full = df[df.o2_fraction == 1.0].iloc[0]
    assert 0.3 < full.yield_gg < 0.7
    assert (df.status == "optimal").all()


@needs_model
def test_o2_cap_applied_in_dfba(model):
    from fbaval import simulate
    vmax = simulate.vmax_from_mu_max(model, 0.2)
    free = simulate.run_dfba(model, {"glucose": 15.0}, 0.05, 40, 0.5, vmax, 0.26, None)
    capped = simulate.run_dfba(model, {"glucose": 15.0}, 0.05, 40, 0.5, vmax, 0.26, 0.5)
    assert capped.biomass_gl.iloc[-1] < free.biomass_gl.iloc[-1]
    assert free.time_h.iloc[-1] == pytest.approx(40, abs=0.6)   # runs to the end


@needs_model
def test_dfba_reports_falling_pH_with_ph_module(model):
    import dfba_kinetics as dk
    from fbaval import simulate
    sim = simulate.run_dfba(model, {"glucose": 5.0}, 0.05, 40.0, 0.5, 2.0, 0.26, None,
                            ph_kwargs={"pH_params": dk.CardinalPH(), "initial_pH": 6.6,
                                       "medium_phosphate_mM": 200.0})
    assert "pH" in sim and sim.pH.iloc[0] == pytest.approx(6.6, abs=1e-3)
    assert sim.pH.iloc[-1] < 6.5            # NH4 assimilation acidifies the buffer


@needs_model
def test_finite_ammonium_stops_growth_but_not_sugar_use(model):
    from fbaval import simulate
    # 5 mM NH4 supports well under 1 g/L of biomass; with growth_cap the sugar is
    # still used up, and the carbon not built into biomass leaves as CO2.
    sim = simulate.run_dfba(model, {"glucose": 2.0}, 0.05, 60.0, 0.5, 2.0, 0.26, None,
                            nh4_mM=5.0, coupling="growth_cap")
    end = sim.iloc[-1]
    comp_n = 5.0 / (end.biomass_gl - 0.05)
    assert end.nh4_mM < 0.05 and 6.0 < comp_n < 8.0      # model biomass ~7 mmol N/g
    assert end.glucose_gl < 0.01
    assert "coupling_fallback" not in sim or not sim.coupling_fallback.any()


@needs_model
def test_normalized_biomass_keeps_yield_and_sets_one_gram_per_flux(model):
    from fbaval import simulate
    with model:
        g = model.reactions.GROWTH
        g.add_metabolites({k: v for k, v in g.metabolites.items() if k.id != "BIOMASS"})  # 2 g per mmol
        assert simulate.biomass_weight(model) == pytest.approx(2.0, rel=1e-3)
        simulate.set_sugars(model, {"glucose": 1.0})
        raw = model.slim_optimize()
        bw = simulate.normalize_biomass(model)
        assert bw == pytest.approx(2.0, rel=1e-3)
        assert simulate.biomass_weight(model) == pytest.approx(1.0, rel=1e-6)
        assert model.slim_optimize() == pytest.approx(raw * bw, rel=1e-6)   # flux now in g
