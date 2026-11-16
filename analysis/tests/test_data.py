import pytest

from fbaval import data

needs_presens = pytest.mark.skipif(
    not (data.condition_dir("Glc15") / "Processed_Presens_Data.xlsx").exists(),
    reason="online (PreSens) series not found in data/flask")

needs_vinestock = pytest.mark.skipif(
    not (data.condition_dir("Glc15") / "workbook.xlsx").exists(),
    reason="flask workbooks not found in data/flask")


def test_otr_conversion():
    assert data.VINESTOCK_OTR_MAX_MMOL_L_H == pytest.approx(46.875)


@needs_vinestock
def test_glucose_hplc_matches_known_values():
    df = data.load_vinestock_hplc("Glc15")
    t0 = df[(df.time_h == 0) & (df.compound == "glucose")].sort_values("flask")
    assert list(t0.conc_gl.round(2)) == [17.13, 16.5, 13.14]
    gly = df[df.compound == "glycerol"]
    assert gly.conc_gl.max() == pytest.approx(0.657, abs=1e-3)
    assert set(gly.flask) == {1, 2, 3}          # duplicated 'Flask A1' headers fixed


@needs_vinestock
def test_final_biomass():
    fin = data.load_vinestock_final("Glc15")
    assert list(fin.biomass_gl.round(2)) == [7.28, 7.33, 7.72]
    assert fin.yield_gg.between(0.4, 0.6).all()


@needs_vinestock
def test_fructose_has_no_glycerol():
    df = data.load_vinestock_hplc("Fru15")
    assert df[df.compound == "glycerol"].conc_gl.max() == 0


def test_four_conditions_with_15_gl_sugar():
    assert set(data.PRIMARY_VINESTOCK) == set(data.VINESTOCK)
    for spec in data.VINESTOCK.values():
        assert spec["glucose_gl"] + spec["fructose_gl"] == 15.0


@needs_presens
@pytest.mark.parametrize("cond", list(data.VINESTOCK))
def test_presens_time_series_long_format(cond):
    df = data.load_vinestock_presens(cond)
    assert list(df.columns) == ["time_h", "flask", "biomass_gl", "ph", "do_pct_as", "our_g_l_h"]
    assert set(df.flask) == {1, 2, 3}
    assert 40 < df.time_h.max() < 70
    assert df.ph.dropna().between(5.8, 6.8).all()


@needs_presens
@pytest.mark.parametrize("cond, first_reading_h", [("Glc15", 0.138), ("Glc10Fru5", 1.009), ("Glc5Fru10", -0.006), ("Fru15", 0.323)])
def test_presens_time_is_measured_from_the_first_hplc_sample(cond, first_reading_h):
    # HPLC and online times must share one zero; the readers started after TP0 in three of the four runs.
    assert data.load_vinestock_presens(cond).time_h.min() == pytest.approx(first_reading_h, abs=0.002)


@needs_presens
@needs_vinestock
@pytest.mark.parametrize("cond", list(data.VINESTOCK))
def test_presens_biomass_ends_at_dry_weight(cond):
    # The reflectance signal is calibrated to the end-point dry weight. Single
    # flasks differ by up to 7 % (Fru15), so compare condition means.
    df = data.load_vinestock_presens(cond)
    last = df.sort_values("time_h").groupby("flask").biomass_gl.last().mean()
    fin = data.load_vinestock_final(cond).biomass_gl.mean()
    assert abs(last - fin) / fin < 0.05
