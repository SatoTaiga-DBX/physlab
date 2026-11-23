"""二体問題: 力 = −∇U、初期条件（ケプラー軌道）、保存量。"""

import numpy as np
import pytest

from physlab.domain.mechanics.two_body import (
    TwoBodyParams,
    accel,
    center_of_mass,
    energy,
    initial_state_cm,
    potential,
    total_momentum,
)

P = TwoBodyParams()


def test_force_is_minus_gradient_of_potential():
    """m_a a_a = −∂U/∂r_a を、U の数値微分（独立な計算）で確かめる。"""
    q = np.array([0.3, -0.2, -0.9, 0.7])
    mass = np.array([P.m1, P.m1, P.m2, P.m2])
    h = 1e-6
    numeric = np.empty(4)
    for i in range(4):
        plus, minus = q.copy(), q.copy()
        plus[i] += h
        minus[i] -= h
        numeric[i] = -(potential(plus, P) - potential(minus, P)) / (2 * h)
    assert mass * accel(P)(0.0, q) == pytest.approx(numeric, rel=1e-6)


def test_newtons_third_law_total_force_is_zero():
    q = np.array([1.0, 2.0, -0.5, 0.3])
    a = accel(P)(0.0, q)
    force_total = P.m1 * a[0:2] + P.m2 * a[2:4]
    assert force_total == pytest.approx(0.0, abs=1e-14)


def test_initial_state_is_a_bound_kepler_orbit_in_the_cm_frame():
    q0, v0 = initial_state_cm(P)
    assert center_of_mass(q0, P) == pytest.approx(0.0, abs=1e-15)
    assert total_momentum(v0, P) == pytest.approx(0.0, abs=1e-15)
    assert np.linalg.norm(q0[0:2] - q0[2:4]) == pytest.approx(P.periapsis)
    v_rel = v0[0:2] - v0[2:4]
    assert v_rel @ v_rel == pytest.approx(6.0)  # G M (1+e) / r_p = 4 × 1.5
    # E = −G m1 m2 / (2a) （a は相対軌道の長半径）
    assert energy(q0, v0, P) == pytest.approx(-P.G * P.m1 * P.m2 / (2 * P.semi_major_axis))


def test_kepler_third_law_period():
    assert P.period == pytest.approx(2 * np.pi * np.sqrt(2.0))
    assert P.semi_major_axis == pytest.approx(2.0)
