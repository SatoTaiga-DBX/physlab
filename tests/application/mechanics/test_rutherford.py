"""ラザフォード散乱: 数える実験と公式、電荷を広げた場合の違い。"""

import numpy as np
import pytest

from physlab.application.mechanics.rutherford import (
    HEAD_ON,
    ball_cross_section,
    ball_curve,
    energy_transfer,
    max_ball_angle,
    point_charge_cross_section,
)
from physlab.application.mechanics.scattering import solid_angles
from physlab.domain.mechanics.scattering import (
    ball_graze_angle,
    ball_small_angle,
    rutherford_angle,
    rutherford_impact_parameter,
)


def _formula_per_bin(edges: np.ndarray) -> np.ndarray:
    """区間に入る ρ の円環の面積を立体角で割った理論値。"""
    rho = rutherford_impact_parameter(np.minimum(edges, np.pi - 1e-12), HEAD_ON)
    return -np.pi * np.diff(rho**2) / solid_angles(edges)


@pytest.fixture(scope="module")
def counted():
    return point_charge_cross_section()


def test_counting_reproduces_the_inverse_fourth_power_of_sine(counted) -> None:
    theory = _formula_per_bin(counted.edges)
    assert np.max(np.abs(counted.measured - theory) / counted.error) < 4.0
    # 20° 付近と 175° 付近で 2 桁以上違う
    assert counted.measured[0] / counted.measured[-1] > 100


def test_head_on_distance_is_one() -> None:
    assert HEAD_ON == pytest.approx(1.0)


def test_a_ball_smaller_than_the_closest_approach_looks_like_a_point() -> None:
    curve = ball_curve(0.8)
    assert curve.chi == pytest.approx(rutherford_angle(curve.rho, HEAD_ON), abs=1e-7)


def test_a_larger_ball_agrees_only_below_the_graze_angle() -> None:
    curve = ball_curve(1.4)
    point = rutherford_angle(curve.rho, HEAD_ON)
    graze = ball_graze_angle(1.4)
    outside = point < graze - 1e-3
    assert curve.chi[outside] == pytest.approx(point[outside], abs=1e-7)
    assert np.max(np.abs(curve.chi - point)[point > graze]) > 0.1


def test_a_wide_ball_has_no_large_angle_scattering() -> None:
    assert np.degrees(max_ball_angle(3.0)) == pytest.approx(24.27, abs=0.05)
    counted = ball_cross_section(3.0)
    assert counted.counts[counted.centers > np.radians(35)].sum() == 0
    assert counted.counts.sum() > 0


def test_very_wide_ball_approaches_the_projected_charge_formula() -> None:
    rho = np.linspace(0.01, 1.2, 20_000) * 1000.0
    estimate = ball_small_angle(rho, HEAD_ON, 1000.0).max()
    assert max_ball_angle(1000.0) == pytest.approx(estimate, rel=2e-3)


def test_energy_transfer_goes_as_inverse_square() -> None:
    result = energy_transfer()
    assert np.max(np.abs(result.measured - result.expected) / result.error) < 4.0
    assert result.maximum == pytest.approx(2 * (2 / 3) ** 2 / 2)
    slope = np.polyfit(np.log(result.centers), np.log(result.measured), 1)[0]
    assert slope == pytest.approx(-2.0, abs=0.05)
