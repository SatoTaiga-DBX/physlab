"""対称性と保存量の実験: 外場で対称性を壊す、積分法の比較、回転する棒の上の質点。

設計書: docs/design/006-conservation-laws.md ／ 教科書: 第 4 章「保存法則」
返すのは dataclass（配列と派生量）だけ。これを描画するのは infrastructure の仕事。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid

from physlab.domain.mechanics.cluster import (
    FIELDS,
    ClusterParams,
    accel,
    angular_momentum,
    center_of_mass,
    energy,
    momentum,
    predicted_rates,
)
from physlab.domain.mechanics.rotating_rod import (
    energy_function,
    kinetic_energy,
    motor_power,
    radius,
)
from physlab.numerics.integrators import Array, integrate, velocity_verlet

# 3 質点の既定の初期状態（全運動量も角運動量も 0 でないようにとる）
DEFAULT_Q0 = (0.0, 0.0, 1.2, 0.0, 0.4, 0.9)
DEFAULT_V0 = (0.3, -0.2, 0.1, 0.8, -0.9, 0.1)


@dataclass(frozen=True)
class ClusterRun:
    """1 つの外場での時間発展と、対称性の破れから予測した保存量の変化。"""

    params: ClusterParams
    t: Array
    q: Array  # (len(t), 2N)
    v: Array
    energy: Array
    momentum: Array  # (len(t), 2)
    angular: Array
    predicted_energy: Array  # E(0) + ∫ (∂U_ext/∂t) dt
    predicted_momentum: Array  # P(0) + ∫ ΣF_ext dt
    predicted_angular: Array  # L(0) + ∫ Σ(r × F_ext) dt


def run_cluster(
    params: ClusterParams,
    q0: Array = DEFAULT_Q0,
    v0: Array = DEFAULT_V0,
    duration: float = 10.0,
    dt: float = 1e-3,
    method: str = "verlet",
) -> ClusterRun:
    """``method``（verlet / rk4）で解き、保存量とその予測値を求める。"""
    t = np.arange(0.0, duration + dt / 2, dt)
    q0, v0 = np.asarray(q0, dtype=float), np.asarray(v0, dtype=float)
    acc = accel(params)
    if method == "verlet":
        q, v = velocity_verlet(acc, q0, v0, t)
    else:
        n = q0.size

        def f(time: float, y: Array) -> Array:
            return np.concatenate([y[n:], acc(time, y[:n])])

        y = integrate(f, np.concatenate([q0, v0]), t, method=method)
        q, v = y[:, :n], y[:, n:]

    e = energy(t, q, v, params)
    p = momentum(v, params)
    ang = angular_momentum(q, v, params)
    d_e, d_p, d_l = predicted_rates(t, q, params)
    return ClusterRun(
        params=params,
        t=t,
        q=q,
        v=v,
        energy=e,
        momentum=p,
        angular=ang,
        predicted_energy=e[0] + cumulative_trapezoid(d_e, t, initial=0.0),
        predicted_momentum=p[0] + cumulative_trapezoid(d_p, t, axis=0, initial=0.0),
        predicted_angular=ang[0] + cumulative_trapezoid(d_l, t, initial=0.0),
    )


def symmetry_experiment(duration: float = 10.0, dt: float = 1e-3) -> dict[str, ClusterRun]:
    """同じ初期状態を、4 通りの外場で解く。"""
    return {
        field: run_cluster(ClusterParams(field=field), duration=duration, dt=dt) for field in FIELDS
    }


def center_of_mass_invariant(run: ClusterRun) -> Array:
    """``M R(t) − P(t) t``。自由な系では一定（ガリレイ変換の対称性に対応する保存量）。"""
    mass = run.params.total_mass
    return mass * center_of_mass(run.q, run.params) - run.momentum * run.t[:, None]


# ---- 積分法と保存量 ----------------------------------------------------------------------


@dataclass(frozen=True)
class IntegratorComparison:
    """自由な系（全部が保存するはず）で、積分法ごとの保存量のずれ。"""

    dt: float
    t: Array
    energy_error: dict[str, Array]  # |E − E0|
    momentum_error: dict[str, Array]  # max(|Px − Px0|, |Py − Py0|)
    angular_error: dict[str, Array]  # |L − L0|


def compare_integrators(duration: float = 50.0, dt: float = 0.01) -> IntegratorComparison:
    errors: tuple[dict, dict, dict] = ({}, {}, {})
    t = np.empty(0)
    for method in ("verlet", "rk4"):
        run = run_cluster(ClusterParams(), duration=duration, dt=dt, method=method)
        t = run.t
        errors[0][method] = np.abs(run.energy - run.energy[0])
        errors[1][method] = np.max(np.abs(run.momentum - run.momentum[0]), axis=1)
        errors[2][method] = np.abs(run.angular - run.angular[0])
    return IntegratorComparison(dt, t, *errors)


# ---- 回転する棒の上の質点 -----------------------------------------------------------------


@dataclass(frozen=True)
class RodRun:
    r0: float
    omega: float
    t: Array
    r: Array
    r_dot: Array
    kinetic: Array  # T（保存しない）
    energy_function: Array  # E = p ṙ − L（保存する）
    motor_work: Array  # モーターがした仕事 ∫ 2 m r ṙ ω² dt
    exact_r: Array


def rotating_rod_run(
    r0: float = 0.5, omega: float = 1.0, duration: float = 1.5, dt: float = 1e-3
) -> RodRun:
    """``r̈ = ω² r`` を RK4 で解く（``ṙ(0) = 0``）。"""
    t = np.arange(0.0, duration + dt / 2, dt)

    def f(_t: float, y: Array) -> Array:
        return np.array([y[1], omega**2 * y[0]])

    y = integrate(f, [r0, 0.0], t)
    r, r_dot = y[:, 0], y[:, 1]
    return RodRun(
        r0=r0,
        omega=omega,
        t=t,
        r=r,
        r_dot=r_dot,
        kinetic=kinetic_energy(r, r_dot, omega),
        energy_function=energy_function(r, r_dot, omega),
        motor_work=cumulative_trapezoid(motor_power(r, r_dot, omega), t, initial=0.0),
        exact_r=radius(t, r0, omega)[0],
    )
