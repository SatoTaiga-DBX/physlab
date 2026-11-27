"""弾性衝突の運動学。"""

import numpy as np
import pytest

from physlab.domain.mechanics.collision import (
    energy_transfer_fraction,
    final_speeds,
    final_velocities,
    lab_angles,
    log_energy_decrement,
    max_deflection,
    momentum_circle,
)

CHI = np.linspace(0.05, np.pi - 0.05, 9)


def _after(m1: float, m2: float, chi: np.ndarray, speed: float = 2.0):
    direction = np.column_stack([np.cos(chi), np.sin(chi)])
    return final_velocities(m1, m2, [speed, 0.0], [0.0, 0.0], direction)


@pytest.mark.parametrize(("m1", "m2"), [(1.0, 3.0), (2.0, 2.0), (4.0, 1.0)])
def test_momentum_and_energy_are_conserved(m1: float, m2: float) -> None:
    u1, u2 = _after(m1, m2, CHI)
    assert m1 * u1 + m2 * u2 == pytest.approx(np.tile([2.0 * m1, 0.0], (CHI.size, 1)))
    kinetic = 0.5 * m1 * np.sum(u1**2, axis=1) + 0.5 * m2 * np.sum(u2**2, axis=1)
    assert kinetic == pytest.approx(0.5 * m1 * 4.0)


@pytest.mark.parametrize(("m1", "m2"), [(1.0, 3.0), (2.0, 2.0), (4.0, 1.0)])
def test_closed_formulas_match_the_vector_construction(m1: float, m2: float) -> None:
    u1, u2 = _after(m1, m2, CHI)
    theta1, theta2 = lab_angles(CHI, m1, m2)
    assert theta1 == pytest.approx(np.arctan2(u1[:, 1], u1[:, 0]))
    assert theta2 == pytest.approx(np.abs(np.arctan2(u2[:, 1], u2[:, 0])))
    s1, s2 = final_speeds(CHI, m1, m2, 2.0)
    assert s1 == pytest.approx(np.linalg.norm(u1, axis=1))
    assert s2 == pytest.approx(np.linalg.norm(u2, axis=1))
    assert energy_transfer_fraction(CHI, m1, m2) == pytest.approx(
        m2 * np.sum(u2**2, axis=1) / (m1 * 4.0)
    )
    center, radius = momentum_circle(m1, m2, 2.0)
    assert np.hypot(m2 * u2[:, 0] - center, m2 * u2[:, 1]) == pytest.approx(radius)


def test_equal_masses_leave_at_right_angles() -> None:
    theta1, theta2 = lab_angles(CHI, 1.0, 1.0)
    assert theta1 + theta2 == pytest.approx(np.pi / 2)
    assert theta1 == pytest.approx(CHI / 2)


def test_heavy_projectile_has_a_maximum_deflection() -> None:
    theta1, _ = lab_angles(np.linspace(0.0, np.pi, 200_001), 4.0, 1.0)
    assert theta1.max() == pytest.approx(max_deflection(4.0, 1.0), abs=1e-8)
    assert max_deflection(4.0, 1.0) == pytest.approx(np.arcsin(0.25))
    assert max_deflection(1.0, 4.0) == pytest.approx(np.pi)


def test_energy_transfer_is_complete_only_for_equal_masses() -> None:
    assert energy_transfer_fraction(np.pi, 1.0, 1.0) == pytest.approx(1.0)
    assert energy_transfer_fraction(np.pi, 1.0, 12.0) == pytest.approx(48.0 / 169.0)


def test_log_energy_decrement() -> None:
    assert log_energy_decrement(1.0) == 1.0
    # 重い標的では ξ ≈ 2 / (A + 2/3)
    assert log_energy_decrement(100.0) == pytest.approx(2.0 / (100.0 + 2.0 / 3.0), rel=1e-4)
    # cos χ を一様に選んだ平均と一致する（A = 12）
    cos_chi = np.linspace(-1.0, 1.0, 200_001)
    ratio = (1.0 + 12.0**2 + 2.0 * 12.0 * cos_chi) / 13.0**2
    average = np.trapezoid(-np.log(ratio), cos_chi) / 2.0
    assert log_energy_decrement(12.0) == pytest.approx(average, rel=1e-8)
