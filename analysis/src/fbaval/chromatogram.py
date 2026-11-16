"""Area of the refractive-index peak at 16.1 min, which elutes on the tail of the fructose peak.

The peak table of the instrument software integrates it inconsistently on that tail. Here a chromatogram of the same
fermentation before any product formed is scaled and shifted to match the tail ahead of the peak, then subtracted.
The remainder is integrated above a straight baseline across the whole peak.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

TAIL_WINDOW = (15.20, 15.60)      # min, only the fructose tail
PEAK_WINDOW = (15.62, 16.92)      # min, the whole peak
LEFT_BASELINE = (15.56, 15.62)
RIGHT_BASELINE = (16.92, 17.05)
ZERO_WINDOW = (17.05, 17.25)      # min, flat signal that sets each trace's zero


def _zeroed(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    return y - np.median(y[(t > ZERO_WINDOW[0]) & (t < ZERO_WINDOW[1])])


def peak_16p1_area(t: np.ndarray, y: np.ndarray, reference: tuple[np.ndarray, np.ndarray] | None = None) -> float:
    """Area in signal x s. ``reference`` is the (time, signal) of a sample without product, given when fructose is present."""
    signal = _zeroed(t, y)
    if reference is not None:
        tr, yr = reference
        yr = _zeroed(tr, yr)
        w = (t >= TAIL_WINDOW[0]) & (t <= TAIL_WINDOW[1])
        fit = least_squares(lambda p: p[0] * np.interp(t[w] - p[1], tr, yr) - signal[w], [1.0, 0.0], bounds=([0.0, -0.3], [5.0, 0.3]))
        signal = signal - fit.x[0] * np.interp(t - fit.x[1], tr, yr)
    left = (t >= LEFT_BASELINE[0]) & (t < LEFT_BASELINE[1])
    right = (t > RIGHT_BASELINE[0]) & (t <= RIGHT_BASELINE[1])
    mid = (t >= PEAK_WINDOW[0]) & (t <= PEAK_WINDOW[1])
    x0, x1, y0, y1 = t[left].mean(), t[right].mean(), signal[left].mean(), signal[right].mean()
    return float(np.trapezoid(signal[mid] - (y0 + (y1 - y0) * (t[mid] - x0) / (x1 - x0)), t[mid] * 60.0))
