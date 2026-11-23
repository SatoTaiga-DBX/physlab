"""二重振り子の物理: ラグランジアン・運動方程式・エネルギー・微小振動の基準振動数。

設計書: docs/design/005-analytical-mechanics-ch1.md ／ 教科書: 第 3 章「運動方程式」§一般座標

一般座標は 2 つの角度 θ1（支点から第 1 の錘）と θ2（第 1 の錘から第 2 の錘）。
どちらも鉛直下向きから測る。
Δ = θ1 − θ2 として::

    T = ½ (m1+m2) l1² ω1² + ½ m2 l2² ω2² + m2 l1 l2 ω1 ω2 cos Δ
    U = −(m1+m2) g l1 cos θ1 − m2 g l2 cos θ2

オイラー＝ラグランジュ方程式は質量行列 M(θ) を使って ``M(θ) θ̈ = F(θ, ω)`` と書ける。
状態は ``y = [θ1, θ2, ω1, ω2]``。
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


@dataclass(frozen=True)
class DoublePendulumParams:
    g: float = 9.81  # 重力加速度 [m/s²]
    m1: float = 1.0  # 第 1 の錘の質量 [kg]
    m2: float = 1.0  # 第 2 の錘の質量 [kg]
    l1: float = 1.0  # 第 1 の棒の長さ [m]
    l2: float = 1.0  # 第 2 の棒の長さ [m]
    theta1: float = np.deg2rad(90.0)  # 初期角 [rad]
    theta2: float = 0.0
    omega1: float = 0.0  # 初期角速度 [rad/s]
    omega2: float = 0.0

    def initial_state(self) -> Array:
        return np.array([self.theta1, self.theta2, self.omega1, self.omega2])


def mass_matrix(theta1: Array, theta2: Array, p: DoublePendulumParams) -> Array:
    """質量行列 ``M(θ)``（``T = ½ ω·M ω``）。スカラーなら形状 ``(2, 2)``。"""
    off = p.m2 * p.l1 * p.l2 * np.cos(theta1 - theta2)
    m = np.empty(np.shape(theta1) + (2, 2))
    m[..., 0, 0] = (p.m1 + p.m2) * p.l1**2
    m[..., 0, 1] = m[..., 1, 0] = off
    m[..., 1, 1] = p.m2 * p.l2**2
    return m


def kinetic_energy(
    theta1: Array, theta2: Array, omega1: Array, omega2: Array, p: DoublePendulumParams
) -> Array:
    return (
        0.5 * (p.m1 + p.m2) * p.l1**2 * omega1**2
        + 0.5 * p.m2 * p.l2**2 * omega2**2
        + p.m2 * p.l1 * p.l2 * omega1 * omega2 * np.cos(theta1 - theta2)
    )


def potential_energy(theta1: Array, theta2: Array, p: DoublePendulumParams) -> Array:
    return -(p.m1 + p.m2) * p.g * p.l1 * np.cos(theta1) - p.m2 * p.g * p.l2 * np.cos(theta2)


def lagrangian(
    theta1: Array, theta2: Array, omega1: Array, omega2: Array, p: DoublePendulumParams
) -> Array:
    """``L = T − U``。"""
    return kinetic_energy(theta1, theta2, omega1, omega2, p) - potential_energy(theta1, theta2, p)


def energy(
    theta1: Array, theta2: Array, omega1: Array, omega2: Array, p: DoublePendulumParams
) -> Array:
    """力学的エネルギー ``E = T + U``（保存される）。"""
    return kinetic_energy(theta1, theta2, omega1, omega2, p) + potential_energy(theta1, theta2, p)


def deriv(p: DoublePendulumParams) -> Callable[[float, Array], Array]:
    """``dy/dt = f(t, y)``。``y = [θ1, θ2, ω1, ω2]``。"""

    def f(_t: float, y: Array) -> Array:
        th1, th2, w1, w2 = y
        delta = th1 - th2
        force = np.array(
            [
                -p.m2 * p.l1 * p.l2 * w2**2 * np.sin(delta)
                - (p.m1 + p.m2) * p.g * p.l1 * np.sin(th1),
                p.m2 * p.l1 * p.l2 * w1**2 * np.sin(delta) - p.m2 * p.g * p.l2 * np.sin(th2),
            ]
        )
        accel = np.linalg.solve(mass_matrix(th1, th2, p), force)
        return np.array([w1, w2, accel[0], accel[1]])

    return f


def normal_modes(p: DoublePendulumParams) -> tuple[Array, Array]:
    """微小振動の基準振動 ``(ω_i, 振幅比のベクトル)``。``det(K − ω² M0) = 0`` の根。

    返す ``ω`` は昇順。``vectors[:, i]`` が ``i`` 番目の基準振動の ``(θ1, θ2)`` の
    振幅比（正規化済み）。
    """
    m0 = mass_matrix(np.float64(0.0), np.float64(0.0), p)
    k = np.diag([(p.m1 + p.m2) * p.g * p.l1, p.m2 * p.g * p.l2])
    eigenvalues, vectors = np.linalg.eig(np.linalg.solve(m0, k))
    order = np.argsort(eigenvalues.real)
    omega = np.sqrt(eigenvalues.real[order])
    vectors = vectors.real[:, order]
    return omega, vectors / np.linalg.norm(vectors, axis=0)


def positions(
    theta1: Array, theta2: Array, p: DoublePendulumParams
) -> tuple[Array, Array, Array, Array]:
    """2 つの錘の位置 ``(x1, y1, x2, y2)``。支点が原点、``y`` は鉛直上向き。"""
    x1 = p.l1 * np.sin(theta1)
    y1 = -p.l1 * np.cos(theta1)
    return x1, y1, x1 + p.l2 * np.sin(theta2), y1 - p.l2 * np.cos(theta2)
