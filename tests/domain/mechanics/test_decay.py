"""崩壊の運動学: 保存則と、角度分布の性質。"""

import numpy as np
import pytest

from physlab.domain.mechanics.decay import (
    DecayParams,
    angle_density,
    kinetic_energy_range,
    lab_angle,
    lab_speeds,
    max_energy_in_many_body,
    max_lab_angle,
)


def test_rest_frame_momenta_share_the_released_energy() -> None:
    p = DecayParams(m1=1.0, m2=3.0, energy=24.0)
    v1, v2 = p.rest_speeds
    assert p.momentum == pytest.approx(6.0)
    assert (v1, v2) == pytest.approx((6.0, 2.0))
    assert 0.5 * p.m1 * v1**2 + 0.5 * p.m2 * v2**2 == pytest.approx(p.energy)


def test_lab_speeds_solve_the_velocity_triangle() -> None:
    v0, theta0 = 2.0, np.linspace(0.1, 3.0, 7)
    for big_v in (1.0, 3.0):
        k = big_v / v0
        theta = lab_angle(theta0, k)
        speed = np.hypot(v0 * np.cos(theta0) + big_v, v0 * np.sin(theta0))
        fast, slow = lab_speeds(theta, big_v, v0)
        matches = np.minimum(np.abs(fast - speed), np.nan_to_num(np.abs(slow - speed), nan=1e9))
        assert np.max(matches) < 1e-9


def test_fast_parent_confines_the_fragments_to_a_cone() -> None:
    k = 2.0
    theta = lab_angle(np.linspace(0.0, np.pi, 100_001), k)
    assert theta.max() == pytest.approx(max_lab_angle(k), abs=1e-8)
    assert max_lab_angle(k) == pytest.approx(np.pi / 6)
    assert max_lab_angle(0.5) == pytest.approx(np.pi)


def test_density_of_a_parent_at_rest_is_isotropic() -> None:
    theta = np.linspace(0.0, np.pi, 11)
    assert angle_density(theta, 0.0) == pytest.approx(0.5 * np.sin(theta), abs=1e-14)


@pytest.mark.parametrize("k", [0.0, 0.5, 0.9, 1.0, 1.5, 3.0])
def test_density_integrates_to_one(k: float) -> None:
    # 上限の角の 1/√ の特異性を、置換で打ち消してから台形公式で積分する
    s = np.linspace(0.0, 1.0, 400_001)[1:]
    top = max_lab_angle(k)
    theta = top - top * s**2
    value = np.trapezoid(angle_density(theta, k) * 2.0 * top * s, s)
    assert value == pytest.approx(1.0, abs=2e-4)


def test_energy_range_and_many_body_maximum() -> None:
    assert kinetic_energy_range(2.0, 1.0, 3.0) == pytest.approx((4.0, 16.0))
    assert max_energy_in_many_body(1.0, 4.0, 10.0) == pytest.approx(7.5)
