import numpy as np
import pytest

from fbaval import ph


def test_rhm_flask_ionic_strength_is_about_0p4_molar():
    # 0.2 M KH2PO4 + 76 mM NH4Cl dominate; Na+ comes from the charge balance.
    assert 0.30 < ph.ionic_strength(6.6) < 0.45


def test_no_growth_keeps_initial_pH():
    t = np.linspace(0, 50, 6)
    out = ph.predict_ph(np.zeros_like(t), {}, 6.6, ph.VERSIONS["v3"], ph.Composition(7.0, 0.4, 0.08))
    assert np.allclose(out, 6.6, atol=1e-3)


def test_ammonium_makes_pH_fall_monotonically_with_biomass():
    x = np.linspace(0, 8, 9)
    out = ph.predict_ph(x, {}, 6.6, ph.VERSIONS["v1"], ph.Composition(7.0, 0.4, 0.08))
    assert np.all(np.diff(out) < 0) and 6.6 - out[-1] > 0.3


def test_v0_ignores_assimilation():
    x = np.linspace(0, 8, 9)
    out = ph.predict_ph(x, {}, 6.6, ph.VERSIONS["v0"], ph.Composition(7.0, 0.4, 0.08))
    assert np.allclose(out, 6.6, atol=1e-3)


def test_measured_acetate_lowers_pH():
    x = np.zeros(3)
    acids = {"acetate": np.array([0.0, 5.0, 10.0])}
    out = ph.predict_ph(x, acids, 6.6, ph.VERSIONS["v0"], ph.Composition(7.0, 0.4, 0.08))
    assert out[0] == pytest.approx(6.6, abs=1e-3) and out[2] < out[1] < out[0]


def test_flask_medium_ionic_strength():
    from fbaval import ph
    assert ph.PHOSPHATE_MM == pytest.approx(202.81, abs=0.05)   # 46 mL of RHM with 30 g/L KH2PO4 in 50 mL
    assert ph.NH4CL_MM == pytest.approx(75.68, abs=0.05)
    assert ph.ionic_strength(6.6) == pytest.approx(0.446, abs=0.002)
