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
