"""二重振り子: 運動方程式がラグランジアンから正しく出ていること、基準振動、極限。"""

import numpy as np
import pytest

from physlab.domain.mechanics.double_pendulum import (
    DoublePendulumParams,
    deriv,
    energy,
    lagrangian,
    mass_matrix,
    normal_modes,
    positions,
)

P = DoublePendulumParams(g=9.81, m1=1.3, m2=0.7, l1=1.1, l2=0.9)


def test_equilibrium_is_stationary_and_energy_is_minus_mgl():
    f = deriv(P)
    assert np.allclose(f(0.0, np.zeros(4)), 0.0)
    e = energy(0.0, 0.0, 0.0, 0.0, P)
    assert e == pytest.approx(-(P.m1 + P.m2) * P.g * P.l1 - P.m2 * P.g * P.l2)


def test_mass_matrix_is_symmetric_positive_definite():
    for th1, th2 in [(0.0, 0.0), (1.0, -2.0), (3.0, 0.1)]:
        m = mass_matrix(np.float64(th1), np.float64(th2), P)
        assert np.allclose(m, m.T)
        assert np.all(np.linalg.eigvalsh(m) > 0)


def test_equations_of_motion_satisfy_euler_lagrange_from_the_lagrangian():
    """d/dt(∂L/∂ω) − ∂L/∂θ = 0 を、L の数値微分だけで確かめる（deriv の式とは独立）。"""
    state = np.array([0.7, -0.4, 1.3, -0.8])
    accel = deriv(P)(0.0, state)[2:]

    def L(y: np.ndarray) -> float:
        return float(lagrangian(y[0], y[1], y[2], y[3], P))

    def grad(f, y: np.ndarray, idx: slice, h: float) -> np.ndarray:
        out = []
        for i in range(idx.start, idx.stop):
            plus, minus = y.copy(), y.copy()
            plus[i] += h
            minus[i] -= h
            out.append((f(plus) - f(minus)) / (2 * h))
        return np.array(out)

    def momentum(y: np.ndarray) -> np.ndarray:  # p = ∂L/∂ω（L は ω の 2 次式なので中心差分が厳密）
        return grad(L, y, slice(2, 4), 1e-3)

    h = 1e-5
    dp_dt = np.zeros(2)
    for i in range(4):  # d/dt p = Σ_i (∂p/∂y_i) ẏ_i, ẏ = (ω, θ̈)
        plus, minus = state.copy(), state.copy()
        plus[i] += h
        minus[i] -= h
        y_dot_i = state[2 + i] if i < 2 else accel[i - 2]
        dp_dt += (momentum(plus) - momentum(minus)) / (2 * h) * y_dot_i
    dL_dtheta = grad(L, state, slice(0, 2), 1e-6)
    assert np.allclose(dp_dt - dL_dtheta, 0.0, atol=1e-5)


def test_normal_modes_of_equal_masses_and_lengths():
    p = DoublePendulumParams(g=9.81, m1=1.0, m2=1.0, l1=0.8, l2=0.8)
    omega, vectors = normal_modes(p)
    w2 = np.array([2 - np.sqrt(2), 2 + np.sqrt(2)]) * p.g / p.l1
    assert omega**2 == pytest.approx(w2, rel=1e-12)
    # 低い方は同位相（θ2/θ1 = √2）、高い方は逆位相（−√2）
    assert vectors[1, 0] / vectors[0, 0] == pytest.approx(np.sqrt(2), rel=1e-12)
    assert vectors[1, 1] / vectors[0, 1] == pytest.approx(-np.sqrt(2), rel=1e-12)


def test_light_second_bob_decouples_into_two_single_pendulums():
    p = DoublePendulumParams(g=9.81, m1=1.0, m2=1e-9, l1=1.2, l2=0.7)
    omega, _ = normal_modes(p)
    expected = np.sort([np.sqrt(p.g / p.l1), np.sqrt(p.g / p.l2)])
    assert omega == pytest.approx(expected, rel=1e-6)


def test_positions_of_hanging_pendulum():
    x1, y1, x2, y2 = positions(0.0, 0.0, P)
    assert (x1, y1) == pytest.approx((0.0, -P.l1))
    assert (x2, y2) == pytest.approx((0.0, -P.l1 - P.l2))
