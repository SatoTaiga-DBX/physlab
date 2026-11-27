"""粒子の崩壊の運動学（保存則だけで決まる部分）。

設計書: docs/design/007-collisions.md ／ 教科書: 第 5 章 §粒子の崩壊

内部エネルギー ``ε`` が解放されて、親が 2 つの破片（質量 ``m1``, ``m2``）に分かれる。
親が静止している系では運動量が打ち消し合うので ``p0² / (2 m1) + p0² / (2 m2) = ε``、
すなわち ``p0 = √(2 μ ε)``（μ は換算質量）。
親が速度 V で動く系での破片の速度は、静止系の速度 ``v0 n`` に V を足したもの。

角度の分布は、速さの比 ``k = V / v0`` だけで決まる。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

Array = NDArray[np.float64]


def reduced_mass(m1: float, m2: float) -> float:
    return m1 * m2 / (m1 + m2)


@dataclass(frozen=True)
class DecayParams:
    m1: float = 1.0
    m2: float = 3.0
    energy: float = 24.0  # 解放されるエネルギー ε（運動エネルギーになる分）

    def __post_init__(self) -> None:
        if min(self.m1, self.m2) <= 0 or self.energy <= 0:
            raise ValueError("質量と解放エネルギーは正の値")

    @property
    def momentum(self) -> float:
        """静止系での 1 つの破片の運動量の大きさ ``p0 = √(2 μ ε)``。"""
        return float(np.sqrt(2.0 * reduced_mass(self.m1, self.m2) * self.energy))

    @property
    def rest_speeds(self) -> tuple[float, float]:
        """静止系での破片の速さ ``(v10, v20)``。"""
        return self.momentum / self.m1, self.momentum / self.m2


def lab_angle(theta0: ArrayLike, k: float) -> Array:
    """静止系での角 θ0（親の進行方向から測る）に対応する、実験室系の角 θ。"""
    theta0 = np.asarray(theta0, dtype=float)
    return np.arctan2(np.sin(theta0), np.cos(theta0) + k)


def max_lab_angle(k: float) -> float:
    """実験室系の角の上限。親が破片より速いか同じ（``k ≥ 1``）と ``arcsin(1/k)``。"""
    return float(np.arcsin(1.0 / k)) if k >= 1.0 else float(np.pi)


def lab_speeds(theta: ArrayLike, v_parent: float, v0: float) -> tuple[Array, Array]:
    """実験室系の向き θ に飛ぶ破片の速さ ``(速い方, 遅い方)``。あり得ない解は NaN。

    速度の三角形 ``v0² = v² + V² − 2 v V cos θ`` を v について解いたもの。
    ``V < v0`` では速い方だけ、``V > v0`` では 2 つとも正になる。
    """
    theta = np.asarray(theta, dtype=float)
    inside = v0**2 - (v_parent * np.sin(theta)) ** 2
    root = np.sqrt(np.where(inside >= 0, inside, np.nan))
    along = v_parent * np.cos(theta)
    fast, slow = along + root, along - root
    return np.where(fast > 0, fast, np.nan), np.where(slow > 0, slow, np.nan)


def angle_density(theta: ArrayLike, k: float) -> Array:
    """静止系で等方的に崩壊するときの、実験室系の角 θ の確率密度（1/rad）。

    ``k = V/v0 < 1`` では θ ごとに速度が 1 つ、``k ≥ 1`` では 2 つの枝を足したもの。
    上限の角 ``arcsin(1/k)`` では、密度が（積分できる形で）無限大になる。
    """
    theta = np.asarray(theta, dtype=float)
    s, c = np.sin(theta), np.cos(theta)
    # 1 − k² sin²θ と 1 + k² cos 2θ を、桁落ちしない形（k = 1 でも安全）で計算する
    inside = (1.0 - k**2) + (k * c) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = ((1.0 - k**2) + 2.0 * (k * c) ** 2) / np.sqrt(np.maximum(inside, 0.0))
        if k < 1.0:
            return 0.5 * s * (2.0 * k * c + ratio)
        allowed = (theta < max_lab_angle(k)) & (c > 0)
        return np.where(allowed, s * ratio, 0.0)


def kinetic_energy_range(mass: float, v_parent: float, v0: float) -> tuple[float, float]:
    """実験室系での破片の運動エネルギーの範囲。等方崩壊ではこの間に一様に分布する。"""
    return 0.5 * mass * (v_parent - v0) ** 2, 0.5 * mass * (v_parent + v0) ** 2


def max_energy_in_many_body(mass: float, parent_mass: float, energy: float) -> float:
    """3 体以上に崩壊するとき、質量 ``mass`` の 1 つが静止系で持ち得る最大の運動エネルギー。

    残りがひとかたまりで反対向きに飛ぶとき（2 体崩壊と同じ形）が最大で ``(M − m) ε / M``。
    """
    return (parent_mass - mass) / parent_mass * energy
