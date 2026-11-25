"""面積速度: 中心力の軌道では、同じ時間に掃く面積が等しい。"""

import numpy as np
import pytest

from physlab.application.mechanics.areal_velocity import central_orbit, equal_time_sectors


@pytest.fixture(scope="module")
def orbit():
    return central_orbit()


def test_angular_momentum_is_conserved_and_energy_stays_bounded(orbit) -> None:
    assert np.max(np.abs(orbit.angular - orbit.angular[0])) < 1e-12
    assert np.max(np.abs(orbit.energy - orbit.energy[0])) < 1e-4


def test_each_verlet_step_sweeps_exactly_the_same_area(orbit) -> None:
    """加速度が r に平行なので、1 ステップの三角形の面積は h |L| / (2m) に丸め誤差で一致する。"""
    h = orbit.t[1] - orbit.t[0]
    expected = h * abs(orbit.angular[0]) / (2 * orbit.params.mass)
    assert orbit.step_area == pytest.approx(expected, rel=1e-9)


def test_equal_time_sectors_have_equal_area_but_very_different_arc_length(orbit) -> None:
    sectors = equal_time_sectors(orbit, 12)
    expected = sectors.duration * abs(orbit.angular[0]) / (2 * orbit.params.mass)
    assert sectors.areas == pytest.approx(expected, rel=1e-9)
    assert sectors.arc_lengths.max() / sectors.arc_lengths.min() > 2.0


def test_orbit_is_a_precessing_rosette_not_a_circle(orbit) -> None:
    radius = np.linalg.norm(orbit.q, axis=1)
    assert radius.max() / radius.min() > 3.0
