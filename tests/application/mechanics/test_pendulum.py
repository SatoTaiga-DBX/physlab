"""単振り子のユースケース（application）: 時間積分した結果を解析解・保存則と照合する。"""

import numpy as np
import pytest

from physlab.application.mechanics.pendulum import (
    METHODS,
    compare_energy_drift,
    compare_with_small_angle,
    measure_period,
    simulate,
)
from physlab.domain.mechanics.pendulum import PendulumParams, energy, period_exact


@pytest.mark.parametrize("theta0_deg", [10, 60, 120, 170])
def test_numeric_period_matches_elliptic_integral(theta0_deg: float) -> None:
    params = PendulumParams(theta0=np.deg2rad(theta0_deg))
    traj = simulate(params, duration=3 * period_exact(params), dt=1e-3, method="rk4")
    assert measure_period(traj) == pytest.approx(period_exact(params), rel=1e-4)


def test_small_angle_limit_trajectory() -> None:
    """振れ角が十分小さければ、非線形解は小角近似の解析解（cos）に一致する。"""
    params = PendulumParams(theta0=np.deg2rad(1.0))
    result = compare_with_small_angle(params, duration=5.0)
    assert np.abs(result.theta - result.theta_linear).max() < 1e-4


def test_large_angle_trajectory_departs_from_small_angle() -> None:
    result = compare_with_small_angle(PendulumParams(theta0=np.deg2rad(120.0)), duration=5.0)
    assert result.period_exact > 1.3 * result.period_linear
    assert np.abs(result.theta - result.theta_linear).max() > 1.0  # rad: 大きくずれる


@pytest.mark.parametrize("method", ["rk4", "verlet"])
def test_energy_conserved(method: str) -> None:
    params = PendulumParams(theta0=np.deg2rad(120.0))
    traj = simulate(params, duration=100.0, dt=0.01, method=method)
    e = energy(traj.theta, traj.omega, params)
    assert np.abs(e / e[0] - 1).max() < 1e-2


def test_energy_comparison_orders_methods_as_documented() -> None:
    """知見 numerical-integration.md: オイラーは増え続け、RK4 が最も正確、ベルレは有界。"""
    params = PendulumParams(theta0=np.deg2rad(120.0))
    result = compare_energy_drift(params, duration=100.0, dt=0.01)
    assert tuple(result.relative_error) == METHODS
    worst = {m: np.abs(err).max() for m, err in result.relative_error.items()}
    assert worst["euler"] > 1.0  # 100% 以上のエネルギー誤差
    assert worst["rk4"] < worst["verlet"] < 1e-3 < worst["euler"]
    assert all(err.shape == result.t.shape for err in result.relative_error.values())


def test_period_ratio_curve_matches_textbook_table_and_series():
    """教科書の表: T/T0 = 1.0174 (30°)、1.3729 (120°)。小振幅では 1 + θ0²/16 に近い。"""
    from physlab.application.mechanics.pendulum import period_ratio_curve

    curve = period_ratio_curve(np.array([10.0, 30.0, 60.0, 120.0, 179.0]))
    assert curve.ratio_exact == pytest.approx([1.0019, 1.0174, 1.0732, 1.3729, 3.9011], abs=5e-5)
    assert curve.ratio_first == pytest.approx([1.0019, 1.0171, 1.0685, 1.2742, 1.6100], abs=5e-5)
    # 2 項まで取ると 60° で誤差 0.04% まで縮まる（教科書の演習の答え）
    assert curve.ratio_second[2] == pytest.approx(1.0728, abs=5e-5)
    assert np.all(np.diff(curve.ratio_exact) > 0)  # 振幅が大きいほど周期は伸びる


def test_phase_portrait_energy_has_its_minimum_at_the_origin_and_separatrix_at_the_top():
    from physlab.application.mechanics.pendulum import phase_portrait

    portrait = phase_portrait(PendulumParams(), n=201)
    i0 = int(np.argmin(np.abs(portrait.omega)))  # ω = 0 の行
    j0 = int(np.argmin(np.abs(portrait.theta)))  # θ = 0 の列
    assert portrait.energy[i0, j0] == pytest.approx(0.0, abs=1e-12)
    assert portrait.energy.min() == pytest.approx(0.0, abs=1e-12)  # E ≥ 0（最小は釣り合いの位置）
    # 真上 (θ=π, ω=0) のエネルギーが分離線の値 2gL
    j_top = int(np.argmin(np.abs(portrait.theta - np.pi)))
    assert portrait.energy[i0, j_top] == pytest.approx(portrait.separatrix_energy, rel=1e-9)
