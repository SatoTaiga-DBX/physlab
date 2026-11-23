"""1 次元の質点のラグランジアン・作用・停留経路。

停留経路がオイラー＝ラグランジュ方程式を満たすこと、離散作用の勾配が正しいこと、
離散作用の停留条件が速度ベルレ法の更新式に一致することを確かめる。
"""

import numpy as np
import pytest

from physlab.domain.mechanics.lagrangian import (
    System1D,
    action,
    discrete_action,
    discrete_action_gradient,
    free_particle,
    harmonic_oscillator,
    sine_mode,
    uniform_gravity,
)

SYSTEMS = [
    pytest.param(free_particle(2.0), 1.0, id="free"),
    pytest.param(uniform_gravity(1.5, 9.81), 2.0, id="gravity"),
    pytest.param(harmonic_oscillator(0.7, 2 * np.pi), 0.3, id="oscillator"),
]


@pytest.mark.parametrize(("system", "tau"), SYSTEMS)
def test_stationary_path_satisfies_endpoints_and_euler_lagrange(system: System1D, tau: float):
    t = np.linspace(0.0, tau, 4001)
    q = system.stationary_path(t, 0.0, 1.0, tau, -0.5)
    assert q[0] == pytest.approx(1.0)
    assert q[-1] == pytest.approx(-0.5)
    # m q̈ + V'(q) = 0 を、差分で近似した q̈ で確かめる（内部の点）
    dt = t[1] - t[0]
    q_ddot = (q[2:] - 2 * q[1:-1] + q[:-2]) / dt**2
    residual = system.mass * q_ddot + system.dpotential(q[1:-1])
    scale = max(np.max(np.abs(system.dpotential(q))), system.mass)
    assert np.max(np.abs(residual)) / scale < 1e-5


def test_oscillator_path_is_undefined_when_sin_omega_tau_is_zero():
    system = harmonic_oscillator(1.0, 2 * np.pi)
    with pytest.raises(ValueError, match="sin"):
        system.stationary_path(np.array([0.0, 0.25, 0.5]), 0.0, 1.0, 0.5, 0.0)


def test_action_of_gravity_parabola_matches_closed_form():
    """S* = −g²τ³/24（m = 1、両端の高さ 0）。"""
    system = uniform_gravity(1.0, 9.81)
    tau = 2.0
    t = np.linspace(0.0, tau, 4001)
    q = system.stationary_path(t, 0.0, 0.0, tau, 0.0)
    assert action(system, t, q) == pytest.approx(-(9.81**2) * tau**3 / 24, rel=1e-6)


def test_sine_mode_vanishes_at_both_ends():
    t = np.linspace(1.0, 3.0, 11)
    eta = sine_mode(t, 1.0, 3.0, n=3)
    assert eta[0] == pytest.approx(0.0, abs=1e-12)
    assert eta[-1] == pytest.approx(0.0, abs=1e-12)


def test_action_curvature_formula_special_values():
    free = free_particle(1.0)
    assert free.action_curvature(2.0, 1) == pytest.approx(0.25 * 2.0 * (np.pi / 2.0) ** 2)
    osc = harmonic_oscillator(1.0, 2 * np.pi)
    assert osc.action_curvature(0.5, 1) == pytest.approx(0.0, abs=1e-12)  # τ = π/ω でちょうど 0
    assert osc.action_curvature(0.3, 1) > 0
    assert osc.action_curvature(0.8, 1) < 0
    assert osc.action_curvature(0.8, 2) > 0


def test_discrete_action_gradient_matches_finite_difference():
    system = harmonic_oscillator(0.9, 3.0)
    t = np.linspace(0.0, 0.7, 13)
    rng = np.random.default_rng(0)
    q = rng.normal(size=t.size)
    grad = discrete_action_gradient(system, t, q)
    assert grad[0] == 0.0 and grad[-1] == 0.0  # 両端は固定
    for k in range(1, t.size - 1):
        h = 1e-5
        plus, minus = q.copy(), q.copy()
        plus[k] += h
        minus[k] -= h
        numeric = (discrete_action(system, t, plus) - discrete_action(system, t, minus)) / (2 * h)
        assert grad[k] == pytest.approx(numeric, rel=1e-6, abs=1e-8)


def test_discrete_stationary_path_of_gravity_is_exact_parabola_at_nodes():
    """一様重力では離散作用の停留点が節点上で放物線に厳密に一致する（V が 1 次式のため）。"""
    system = uniform_gravity(1.0, 9.81)
    t = np.linspace(0.0, 2.0, 41)
    q = system.stationary_path(t, 0.0, 0.0, 2.0, 0.0)
    assert np.max(np.abs(discrete_action_gradient(system, t, q))) < 1e-10


def test_discrete_stationarity_is_the_velocity_verlet_position_update():
    """g_k = 0 を q_{k+1} について解くと q_{k+1} = 2 q_k − q_{k-1} − (Δt²/m) V'(q_k)。"""
    system = harmonic_oscillator(1.3, 2.0)
    t = np.linspace(0.0, 1.0, 21)
    dt = t[1] - t[0]
    rng = np.random.default_rng(1)
    q = rng.normal(size=t.size)
    k = 7
    # 勾配が 0 になる q_{k+1} を、勾配の 1 次式であることを使って解く
    trial = q.copy()
    trial[k + 1] = 0.0
    g0 = discrete_action_gradient(system, t, trial)[k]
    trial[k + 1] = 1.0
    g1 = discrete_action_gradient(system, t, trial)[k]
    q_next = -g0 / (g1 - g0)
    verlet = 2 * q[k] - q[k - 1] - dt**2 / system.mass * system.dpotential(q[k])
    assert q_next == pytest.approx(verlet, rel=1e-10)
