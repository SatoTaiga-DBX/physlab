"""求積: 厳密値の分かる積分。"""

import numpy as np
import pytest

from physlab.numerics.quadrature import (
    gauss_legendre,
    integrate_inverse_sqrt_end,
    integrate_smooth,
)


def test_weights_sum_to_the_interval_length() -> None:
    nodes, weights = gauss_legendre(20)
    assert weights.sum() == pytest.approx(1.0, abs=1e-14)
    assert np.all((nodes > 0) & (nodes < 1))


def test_smooth_integral_is_exact_for_polynomials() -> None:
    value = integrate_smooth(lambda x: 3 * x**2, 1.0, 2.0, n=4)
    assert float(value) == pytest.approx(7.0, abs=1e-13)


def test_inverse_sqrt_singularity_at_the_end_is_removed() -> None:
    value = integrate_inverse_sqrt_end(lambda x: 1.0 / np.sqrt(1.0 - x**2), 0.0, 1.0, n=40)
    assert float(value) == pytest.approx(np.pi / 2, abs=1e-12)


def test_many_intervals_at_once() -> None:
    upper = np.array([1.0, 2.0, 5.0])
    column = upper[:, None]
    value = integrate_inverse_sqrt_end(
        lambda x: 1.0 / np.sqrt(column - x), np.zeros(3), upper, n=40
    )
    assert value == pytest.approx(2.0 * np.sqrt(upper), rel=1e-12)
