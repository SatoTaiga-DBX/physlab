"""一般座標のユースケース: 「位置だけでは運動が決まらない」実験と、二重振り子の時間発展。

設計書: docs/design/005-analytical-mechanics-ch1.md ／ 教科書: 第 3 章「運動方程式」§一般座標
返すのは dataclass（配列と派生量）だけ。これを描画するのは infrastructure の仕事。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from physlab.application.mechanics.pendulum import Trajectory, measure_period, simulate
from physlab.domain.mechanics.double_pendulum import (
    DoublePendulumParams,
    deriv,
    energy,
    normal_modes,
)
from physlab.domain.mechanics.pendulum import PendulumParams
from physlab.numerics.integrators import Array, integrate

# ---- 同じ位置、違う速度 ------------------------------------------------------------------


@dataclass(frozen=True)
class SameCoordinateRuns:
    """同じ初期角 ``θ0`` で、初期角速度だけを変えた単振り子の運動。"""

    params: PendulumParams  # 長さ・重力加速度（θ0, ω0 は個々の実行で変わる）
    t: Array
    theta0: float
    omega0s: Array
    theta: Array  # 形状 (len(omega0s), len(t))
    omega: Array
    turning_angle_measured: Array  # 実測の最大 |θ|
    turning_angle_energy: Array  # エネルギー保存から決まる折り返しの角度 [rad]（回転運動なら NaN）


def same_position_different_velocity(
    theta0: float,
    omega0s: Array,
    duration: float = 6.0,
    dt: float = 1e-3,
    length: float = 1.0,
    g: float = 9.81,
) -> SameCoordinateRuns:
    """``θ(0) = θ0`` は同じで ``θ̇(0)`` だけ違う運動を、RK4 で解く。

    折り返しの角度 ``θ_turn`` は ``cos θ_turn = cos θ0 − L ω0² / (2g)``（エネルギー保存）。
    """
    omega0s = np.asarray(omega0s, dtype=float)
    base = PendulumParams(g=g, length=length, theta0=theta0)
    runs = [
        simulate(PendulumParams(g=g, length=length, theta0=theta0, omega0=w), duration, dt, "rk4")
        for w in omega0s
    ]
    cos_turn = np.cos(theta0) - length * omega0s**2 / (2.0 * g)
    return SameCoordinateRuns(
        params=base,
        t=runs[0].t,
        theta0=theta0,
        omega0s=omega0s,
        theta=np.array([r.theta for r in runs]),
        omega=np.array([r.omega for r in runs]),
        turning_angle_measured=np.array([np.max(np.abs(r.theta)) for r in runs]),
        turning_angle_energy=np.where(
            cos_turn >= -1.0, np.arccos(np.clip(cos_turn, -1, 1)), np.nan
        ),
    )


# ---- 二重振り子 --------------------------------------------------------------------------


@dataclass(frozen=True)
class DoublePendulumRun:
    params: DoublePendulumParams
    t: Array
    theta1: Array
    theta2: Array
    omega1: Array
    omega2: Array
    energy: Array  # 保存されるはずの力学的エネルギー


def simulate_double_pendulum(
    params: DoublePendulumParams, duration: float, dt: float = 1e-3
) -> DoublePendulumRun:
    """二重振り子を RK4 で解く。"""
    t = np.arange(0.0, duration + dt / 2, dt)
    y = integrate(deriv(params), params.initial_state(), t, method="rk4")
    th1, th2, w1, w2 = y.T
    return DoublePendulumRun(params, t, th1, th2, w1, w2, energy(th1, th2, w1, w2, params))


def measure_normal_mode_period(
    params: DoublePendulumParams,
    mode: int,
    amplitude: float = 1e-3,
    duration: float = 12.0,
    dt: float = 1e-3,
) -> tuple[float, float]:
    """``mode``（0 が低い方）の基準振動で始め、周期を測る。``(実測, 理論 2π/ω)`` を返す。

    初期角を基準振動の振幅比にそろえ、速度 0 で放す。振幅が小さいので微小振動の理論と比べられる。
    """
    omega, vectors = normal_modes(params)
    shape = vectors[:, mode] / np.max(np.abs(vectors[:, mode]))
    start = DoublePendulumParams(
        g=params.g,
        m1=params.m1,
        m2=params.m2,
        l1=params.l1,
        l2=params.l2,
        theta1=amplitude * shape[0],
        theta2=amplitude * shape[1],
    )
    run = simulate_double_pendulum(start, duration, dt)
    measured = measure_period(Trajectory(run.t, run.theta1, run.omega1))
    return measured, float(2.0 * np.pi / omega[mode])
