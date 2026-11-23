"""単振り子の物理: 運動方程式・エネルギー・厳密な周期・小角近似解。

設計書: docs/design/002-pendulum.md ／ 教科書: docs/textbook/01-mechanics/01-pendulum.md

運動方程式（質量 m の質点・長さ L の軽い棒・重力加速度 g、摩擦なし）::

    d²θ/dt² = -(g/L) sin θ

小角近似 sin θ ≈ θ では θ(t) = θ0 cos(ω0 t)（ω0 = √(g/L)）となるが、
振れ角が大きいと周期が伸びて近似が破綻する。

最内層（domain）なので、純粋な関数とデータだけを置く。
時間積分の進行は application、描画は infrastructure の仕事。
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.special import ellipk

Array = NDArray[np.float64]


@dataclass(frozen=True)
class PendulumParams:
    g: float = 9.81  # 重力加速度 [m/s²]
    length: float = 1.0  # 振り子の長さ [m]
    theta0: float = np.deg2rad(120.0)  # 初期角（鉛直下向きから）[rad]
    omega0: float = 0.0  # 初期角速度 [rad/s]

    @property
    def small_angle_omega(self) -> float:
        """小角近似の角振動数 ω0 = √(g/L) [rad/s]。"""
        return float(np.sqrt(self.g / self.length))


def deriv(params: PendulumParams) -> Callable[[float, Array], Array]:
    """状態 y = [θ, ω] に対する dy/dt = f(t, y) を返す。"""

    def f(_t: float, y: Array) -> Array:
        return np.array([y[1], -params.g / params.length * np.sin(y[0])])

    return f


def accel(params: PendulumParams) -> Callable[[float, Array], Array]:
    """角加速度 a(t, θ) を返す（速度ベルレ法のような 2 階形式の解法用）。"""

    def a(_t: float, theta: Array) -> Array:
        return -params.g / params.length * np.sin(theta)

    return a


def energy(theta: Array, omega: Array, params: PendulumParams) -> Array:
    """単位質量あたりの力学的エネルギー E/m [J/kg]。最下点を位置エネルギーの基準とする。"""
    kinetic = 0.5 * (params.length * omega) ** 2
    potential = params.g * params.length * (1.0 - np.cos(theta))
    return kinetic + potential


def period_small_angle(params: PendulumParams) -> float:
    """小角近似の周期 T0 = 2π√(L/g) [s]。振幅によらない（等時性）。"""
    return float(2.0 * np.pi / params.small_angle_omega)


def period_exact(params: PendulumParams) -> float:
    """初期角速度 0 のときの厳密な周期 T = (4/ω0) K(sin²(θ0/2)) [s]。

    K は第 1 種完全楕円積分。小角では T → T0 に収束する。
    """
    m = np.sin(params.theta0 / 2) ** 2  # scipy の ellipk は母数 m = k² を取る
    return float(4.0 / params.small_angle_omega * ellipk(m))


def small_angle_solution(t: Array, params: PendulumParams) -> Array:
    """小角近似の解析解 θ(t) = θ0 cos(ω0 t) + (ω_init/ω0) sin(ω0 t)。"""
    w = params.small_angle_omega
    return params.theta0 * np.cos(w * t) + params.omega0 / w * np.sin(w * t)
