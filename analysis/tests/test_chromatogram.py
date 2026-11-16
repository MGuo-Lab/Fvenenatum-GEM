import numpy as np
import pytest

from fbaval.chromatogram import peak_16p1_area

T = np.arange(14.0, 18.0, 1 / 120)     # min, the 0.5 s sampling of the RID traces


def tail(height: float, shift: float = 0.0) -> np.ndarray:
    """Tail of the fructose peak, which elutes just before the window."""
    return height * np.exp(-(T - 14.7 - shift) / 0.18) * (T > 14.7 + shift)


def peak(area: float) -> np.ndarray:
    """A peak near 16.1 min with the given area in signal x s, about as wide as in the fermentation traces (half height 0.38 min)."""
    sigma = 0.16
    return area / (sigma * 60 * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((T - 16.1) / sigma) ** 2)


def test_area_without_fructose_tail():
    assert peak_16p1_area(T, peak(1500.0)) == pytest.approx(1500.0, rel=0.03)


def test_scaled_fructose_tail_is_removed():
    signal = tail(700.0, shift=0.01) + peak(700.0)
    assert peak_16p1_area(T, signal, reference=(T, tail(1000.0))) == pytest.approx(700.0, rel=0.05)


def test_tail_alone_gives_no_area():
    assert peak_16p1_area(T, tail(700.0, shift=0.01), reference=(T, tail(1000.0))) == pytest.approx(0.0, abs=20.0)


def test_constant_offset_does_not_change_the_area():
    assert peak_16p1_area(T, peak(1000.0) + 250.0) == pytest.approx(1000.0, rel=0.03)
