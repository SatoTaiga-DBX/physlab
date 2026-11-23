"""二体問題（万有引力）の物理: ラグランジアン・力・保存量・軌道の初期条件。

設計書: docs/design/005-analytical-mechanics-ch1.md ／ 教科書: 第 3 章 §質点系のラグランジアン

``L = ½ m1 v1² + ½ m2 v2² − U(|r1 − r2|)``、``U = −G m1 m2 / r``。
``U`` が相対位置だけに依存するので、全運動量 ``P = Σ m v`` が保存し、重心は等速直線運動をする。

状態の並べ方: 位置 ``q = [x1, y1, x2, y2]``、速度 ``v`` も同じ並び（最後の軸が 4 成分）。
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


@dataclass(frozen=True)
class TwoBodyParams:
    m1: float = 3.0
    m2: float = 1.0
    G: float = 1.0  # 万有引力定数（無次元の単位系）
    periapsis: float = 1.0  # 近点での相対距離
    eccentricity: float = 0.5  # 離心率（0 で円軌道）

    @property
    def total_mass(self) -> float:
        return self.m1 + self.m2

    @property
    def semi_major_axis(self) -> float:
        return self.periapsis / (1.0 - self.eccentricity)

    @property
    def period(self) -> float:
        """相対軌道の周期 ``2π √(a³ / (G M))``（ケプラーの第 3 法則）。"""
        return float(2.0 * np.pi * np.sqrt(self.semi_major_axis**3 / (self.G * self.total_mass)))


def initial_state_cm(p: TwoBodyParams) -> tuple[Array, Array]:
    """重心系での初期条件 ``(q0, v0)``。近点（x 軸上）から y 方向に出発する。"""
    r_rel = np.array([p.periapsis, 0.0])
    speed = np.sqrt(p.G * p.total_mass * (1.0 + p.eccentricity) / p.periapsis)  # 近点の相対速度
    v_rel = np.array([0.0, speed])
    q0 = np.concatenate([p.m2 / p.total_mass * r_rel, -p.m1 / p.total_mass * r_rel])
    v0 = np.concatenate([p.m2 / p.total_mass * v_rel, -p.m1 / p.total_mass * v_rel])
    return q0, v0


def split(x: Array) -> tuple[Array, Array]:
    """``[x1, y1, x2, y2]`` を 2 つの位置（または速度）ベクトルに分ける。"""
    return x[..., 0:2], x[..., 2:4]


def potential(q: Array, p: TwoBodyParams) -> Array:
    """``U = −G m1 m2 / |r1 − r2|``。"""
    r1, r2 = split(q)
    return -p.G * p.m1 * p.m2 / np.linalg.norm(r1 - r2, axis=-1)


def accel(p: TwoBodyParams) -> Callable[[float, Array], Array]:
    """``a(t, q)``（``m_a r̈_a = −∂U/∂r_a``）を返す。速度ベルレ法に渡す。"""

    def a(_t: float, q: Array) -> Array:
        r1, r2 = split(q)
        d = r2 - r1
        r = np.linalg.norm(d, axis=-1, keepdims=True)
        pull = p.G * d / r**3  # r1 が r2 に引かれる向きの単位質量あたりの量
        return np.concatenate([p.m2 * pull, -p.m1 * pull], axis=-1)

    return a


def kinetic_energy(v: Array, p: TwoBodyParams) -> Array:
    v1, v2 = split(v)
    return 0.5 * p.m1 * np.sum(v1**2, axis=-1) + 0.5 * p.m2 * np.sum(v2**2, axis=-1)


def lagrangian(q: Array, v: Array, p: TwoBodyParams) -> Array:
    """``L = T − U``。"""
    return kinetic_energy(v, p) - potential(q, p)


def energy(q: Array, v: Array, p: TwoBodyParams) -> Array:
    """``E = T + U``（保存される）。"""
    return kinetic_energy(v, p) + potential(q, p)


def total_momentum(v: Array, p: TwoBodyParams) -> Array:
    """``P = m1 v1 + m2 v2``（保存される）。"""
    v1, v2 = split(v)
    return p.m1 * v1 + p.m2 * v2


def center_of_mass(q: Array, p: TwoBodyParams) -> Array:
    """重心 ``R = (m1 r1 + m2 r2) / M``。"""
    r1, r2 = split(q)
    return (p.m1 * r1 + p.m2 * r2) / p.total_mass
