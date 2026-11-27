"""2 粒子の弾性衝突。保存則で決まる部分と、短距離の力で実際に衝突させるための運動方程式。

設計書: docs/design/007-collisions.md ／ 教科書: 第 5 章 §弾性衝突

重心系では、相対速度の大きさは変わらず、向きだけが角 χ 回る。以下の式は静止した標的
（``v2 = 0``）の場合。

状態の並べ方: 位置 ``q = [x1, y1, x2, y2]`` を衝突の数だけ横に並べる（最後の軸が ``4M``）。
複数の衝突を 1 本の状態ベクトルにまとめて、同時に解けるようにしてある。
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray

from physlab.domain.mechanics.decay import reduced_mass
from physlab.domain.mechanics.potentials import CentralPotential

Array = NDArray[np.float64]


def final_velocities(
    m1: float, m2: float, v1: ArrayLike, v2: ArrayLike, direction: ArrayLike
) -> tuple[Array, Array]:
    """衝突後の速度。``direction`` は重心系で粒子 1 が出ていく向き。

    重心の速度はそのまま、相対速度は大きさを保って ``direction`` を向く。
    """
    v1 = np.asarray(v1, dtype=float)
    v2 = np.asarray(v2, dtype=float)
    n = np.asarray(direction, dtype=float)
    n = n / np.linalg.norm(n, axis=-1, keepdims=True)
    total = m1 + m2
    center = (m1 * v1 + m2 * v2) / total
    speed = np.linalg.norm(v1 - v2, axis=-1, keepdims=True)
    return center + m2 / total * speed * n, center - m1 / total * speed * n


def lab_angles(chi: ArrayLike, m1: float, m2: float) -> tuple[Array, Array]:
    """重心系の角 χ（0〜π）に対する実験室系の角 ``(θ1, θ2)``。標的は静止。

    θ1 は入射粒子の偏向、θ2 は跳ね飛ばされた標的の向き（どちらも入射方向から測る大きさ）。
    """
    chi = np.asarray(chi, dtype=float)
    theta1 = np.arctan2(m2 * np.sin(chi), m1 + m2 * np.cos(chi))
    return theta1, 0.5 * (np.pi - chi)


def final_speeds(chi: ArrayLike, m1: float, m2: float, speed: float) -> tuple[Array, Array]:
    """衝突後の速さ ``(v1', v2')``。``speed`` は入射の速さ。"""
    chi = np.asarray(chi, dtype=float)
    v1 = speed * np.sqrt(m1**2 + m2**2 + 2.0 * m1 * m2 * np.cos(chi)) / (m1 + m2)
    v2 = 2.0 * m1 * speed * np.abs(np.sin(0.5 * chi)) / (m1 + m2)
    return v1, v2


def energy_transfer_fraction(chi: ArrayLike, m1: float, m2: float) -> Array:
    """標的に渡るエネルギーの割合 ``E2'/E1 = 4 m1 m2 / (m1+m2)² · sin²(χ/2)``。"""
    chi = np.asarray(chi, dtype=float)
    return 4.0 * m1 * m2 / (m1 + m2) ** 2 * np.sin(0.5 * chi) ** 2


def max_deflection(m1: float, m2: float) -> float:
    """入射粒子の偏向角 θ1 の上限。重い粒子は軽い標的に大きくは曲げられない。"""
    if m1 < m2:
        return float(np.pi)
    return float(np.arcsin(m2 / m1))


def momentum_circle(m1: float, m2: float, speed: float) -> tuple[float, float]:
    """静止標的の衝突後の ``p2'`` が乗る円の ``(中心の x 座標, 半径)``。どちらも ``μ v``。"""
    radius = reduced_mass(m1, m2) * speed
    return radius, radius


def log_energy_decrement(mass_ratio: float) -> float:
    """重心系で等方的に散乱される粒子の、1 回の衝突での ``ln(E/E')`` の平均 ξ。

    ``mass_ratio`` は標的と入射粒子の質量比 A。``α = ((A−1)/(A+1))²`` として
    ``ξ = 1 + α ln α / (1 − α)``（A = 1 では 1）。中性子の減速の目安になる。
    """
    alpha = ((mass_ratio - 1.0) / (mass_ratio + 1.0)) ** 2
    if alpha == 0.0:
        return 1.0
    return float(1.0 + alpha * np.log(alpha) / (1.0 - alpha))


def accel(m1: float, m2: float, potential: CentralPotential) -> Callable[[float, Array], Array]:
    """``U(|r1 − r2|)`` で押し合う 2 粒子の加速度。速度ベルレ法に渡す。"""

    def a(_t: float, q: Array) -> Array:
        state = q.reshape(*q.shape[:-1], -1, 4)
        d = state[..., 0:2] - state[..., 2:4]
        r = np.linalg.norm(d, axis=-1, keepdims=True)
        force_on_1 = -potential.derivative(r) * d / r
        out = np.concatenate([force_on_1 / m1, -force_on_1 / m2], axis=-1)
        return out.reshape(q.shape)

    return a


def energy(q: Array, v: Array, m1: float, m2: float, potential: CentralPotential) -> Array:
    """各衝突の全エネルギー。形は ``(..., M)``。"""
    state = q.reshape(*q.shape[:-1], -1, 4)
    speed = v.reshape(*v.shape[:-1], -1, 4)
    r = np.linalg.norm(state[..., 0:2] - state[..., 2:4], axis=-1)
    kinetic = 0.5 * m1 * np.sum(speed[..., 0:2] ** 2, axis=-1)
    kinetic += 0.5 * m2 * np.sum(speed[..., 2:4] ** 2, axis=-1)
    return kinetic + potential.value(r)
