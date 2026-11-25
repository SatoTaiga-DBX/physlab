"""中心力: 力が中心を向くこと、円軌道の速さ、指数の範囲。"""

import numpy as np
import pytest

from physlab.domain.mechanics.central import CentralParams, accel, circular_speed, energy


def test_acceleration_points_to_the_origin() -> None:
    p = CentralParams(exponent=1.4)
    q = np.array([0.6, -0.8])
    a = accel(p)(0.0, q)
    assert q[0] * a[1] - q[1] * a[0] == pytest.approx(0.0, abs=1e-14)  # r × a = 0
    assert q @ a < 0


def test_acceleration_is_minus_gradient_of_potential() -> None:
    p = CentralParams(alpha=2.0, exponent=0.7, mass=1.5)
    q = np.array([0.9, 0.4])
    h = 1e-6
    zero = np.zeros(2)
    numeric = np.array(
        [-(energy(q + h * e, zero, p) - energy(q - h * e, zero, p)) / (2 * h) for e in np.eye(2)]
    )
    assert p.mass * accel(p)(0.0, q) == pytest.approx(numeric, rel=1e-6)


def test_circular_speed_balances_the_force() -> None:
    p = CentralParams(exponent=1.4)
    radius = 1.3
    v = circular_speed(radius, p)
    a = accel(p)(0.0, np.array([radius, 0.0]))
    assert -a[0] == pytest.approx(v**2 / radius)


@pytest.mark.parametrize("n", [0.0, 2.0, 2.5])
def test_exponent_outside_the_non_collapsing_range_is_rejected(n: float) -> None:
    with pytest.raises(ValueError, match="exponent"):
        CentralParams(exponent=n)
