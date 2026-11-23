"""単振り子のユースケース: 物理(domain)と数値解法(numerics)を組み合わせ、結果データを返す。

設計書: docs/design/002-pendulum.md
返すのは dataclass（時系列と派生量）だけ。これを描画するのは infrastructure の仕事。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from physlab.domain.mechanics.pendulum import (
    PendulumParams,
    accel,
    deriv,
    energy,
    period_exact,
    period_small_angle,
    small_angle_solution,
)
from physlab.numerics.integrators import Array, integrate, velocity_verlet

METHODS = ("euler", "rk4", "verlet")


@dataclass(frozen=True)
class Trajectory:
    """振り子の時系列。"""

    t: Array
    theta: Array
    omega: Array


@dataclass(frozen=True)
class SmallAngleComparison:
    """厳密（非線形）解と小角近似解の比較結果。"""

    params: PendulumParams
    t: Array
    theta: Array  # 非線形方程式を RK4 で解いた角度
    theta_linear: Array  # 小角近似の解析解
    period_exact: float
    period_linear: float


@dataclass(frozen=True)
class EnergyComparison:
    """積分法ごとのエネルギー誤差の推移。"""

    dt: float
    t: Array
    relative_error: dict[str, Array]  # 手法名 -> (E - E0) / E0


def simulate(params: PendulumParams, duration: float, dt: float, method: str = "rk4") -> Trajectory:
    """``method``（euler / rk4 / verlet）で時間 ``duration`` まで積分する。"""
    t = np.arange(0.0, duration + dt / 2, dt)
    if method == "verlet":
        q, v = velocity_verlet(accel(params), [params.theta0], [params.omega0], t)
        return Trajectory(t, q[:, 0], v[:, 0])
    y = integrate(deriv(params), [params.theta0, params.omega0], t, method=method)
    return Trajectory(t, y[:, 0], y[:, 1])


def measure_period(trajectory: Trajectory) -> float:
    """θ(t) のゼロ交差（線形補間）から周期を推定する。θ0 > 0 から静かに放した場合用。"""
    t, theta = trajectory.t, trajectory.theta
    s = np.sign(theta)
    idx = np.nonzero(s[:-1] * s[1:] < 0)[0]  # 符号が変わった区間
    crossings = t[idx] - theta[idx] * (t[idx + 1] - t[idx]) / (theta[idx + 1] - theta[idx])
    return float(2.0 * (crossings[1] - crossings[0]))  # 隣り合う交差の間隔は半周期


def compare_with_small_angle(
    params: PendulumParams, duration: float, dt: float = 1e-3
) -> SmallAngleComparison:
    """非線形解（RK4）と小角近似の解析解を同じ時刻列で計算する。"""
    traj = simulate(params, duration, dt, "rk4")
    return SmallAngleComparison(
        params=params,
        t=traj.t,
        theta=traj.theta,
        theta_linear=small_angle_solution(traj.t, params),
        period_exact=period_exact(params),
        period_linear=period_small_angle(params),
    )


def compare_energy_drift(
    params: PendulumParams, duration: float, dt: float, methods: tuple[str, ...] = METHODS
) -> EnergyComparison:
    """各積分法で同じ条件を解き、エネルギーの相対誤差 (E - E0)/E0 の推移を返す。"""
    e0 = energy(np.array(params.theta0), np.array(params.omega0), params)
    errors: dict[str, Array] = {}
    t = np.empty(0)
    for method in methods:
        traj = simulate(params, duration, dt, method)
        t = traj.t
        errors[method] = (energy(traj.theta, traj.omega, params) - e0) / e0
    return EnergyComparison(dt=dt, t=t, relative_error=errors)
