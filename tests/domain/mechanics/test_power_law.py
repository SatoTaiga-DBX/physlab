"""冪ポテンシャル: 周期の公式（特別な場合と、独立な数値積分）、相似則。"""

import numpy as np
import pytest
from scipy.integrate import quad

from physlab.domain.mechanics.power_law import accel, period_exact, potential, time_scale


def test_harmonic_oscillator_period_does_not_depend_on_amplitude() -> None:
    for a in (0.1, 1.0, 7.0):
        assert period_exact(2.0, a) == pytest.approx(2 * np.pi, rel=1e-12)


def test_linear_potential_period_grows_like_square_root() -> None:
    """U = |x| は一様な力。静かに放して原点に着くまで √(2A)、周期はその 4 倍。"""
    for a in (0.5, 2.0):
        assert period_exact(1.0, a) == pytest.approx(4 * np.sqrt(2 * a), rel=1e-12)


@pytest.mark.parametrize("k", [1.0, 4.0, 6.0])
def test_period_formula_agrees_with_numerical_quadrature(k: float) -> None:
    """T = 4 ∫0^A dx / √(2 (U(A) − U(x)))。端の特異点は x = A (1 − s²) の置き換えで取り除く。"""
    a = 1.7

    def integrand(s: float) -> float:
        x = a * (1.0 - s * s)
        gap = potential(a, k) - potential(x, k)
        return 2.0 * a * s / np.sqrt(2.0 * gap)

    numeric, _ = quad(integrand, 1e-12, 1.0, limit=200)
    assert period_exact(k, a) == pytest.approx(4 * numeric, rel=1e-6)


@pytest.mark.parametrize("k", [1.0, 2.0, 4.0, 6.0])
def test_similarity_time_scale(k: float) -> None:
    assert period_exact(k, 3.0) / period_exact(k, 1.0) == pytest.approx(time_scale(k, 3.0))
    assert time_scale(k, 3.0) == pytest.approx(3.0 ** (1 - k / 2))


def test_force_is_minus_derivative_of_potential() -> None:
    x = np.array([-1.3, -0.2, 0.7, 2.0])
    h = 1e-6
    numeric = -(potential(x + h, 4.0) - potential(x - h, 4.0)) / (2 * h)
    assert accel(4.0)(0.0, x) == pytest.approx(numeric, rel=1e-6)
