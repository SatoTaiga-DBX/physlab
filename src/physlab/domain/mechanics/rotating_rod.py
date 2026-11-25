"""一定の角速度で回る水平な棒に通した質点（時間に依存する拘束の例）。

設計書: docs/design/006-conservation-laws.md ／ 教科書: 第 4 章 §エネルギー

``x = r cos ωt, y = r sin ωt`` より ``L = ½ m (ṙ² + ω² r²)``。``L`` は時刻を陽に含まないので
エネルギー ``E = p ṙ − L = ½ m (ṙ² − ω² r²)`` は保存する。
ところが運動エネルギー ``T = ½ m (ṙ² + ω² r²)`` は保存しない（棒を回し続けるモーターが仕事をする）。
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


def radius(t: Array, r0: float, omega: float) -> tuple[Array, Array]:
    """``r(0) = r0, ṙ(0) = 0`` の解 ``(r, ṙ)``: ``r = r0 cosh ωt``。"""
    return r0 * np.cosh(omega * t), r0 * omega * np.sinh(omega * t)


def kinetic_energy(r: Array, r_dot: Array, omega: float, mass: float = 1.0) -> Array:
    """実際の運動エネルギー ``T = ½ m (ṙ² + ω² r²)``。"""
    return 0.5 * mass * (r_dot**2 + omega**2 * r**2)


def energy_function(r: Array, r_dot: Array, omega: float, mass: float = 1.0) -> Array:
    """ラグランジアンから定まるエネルギー ``E = p ṙ − L = ½ m (ṙ² − ω² r²)``（保存する）。"""
    return 0.5 * mass * (r_dot**2 - omega**2 * r**2)


def motor_power(r: Array, r_dot: Array, omega: float, mass: float = 1.0) -> Array:
    """棒を一定の角速度で回すためにモーターが与える仕事率 ``(dJ/dt) ω = 2 m r ṙ ω²``。"""
    return 2.0 * mass * r * r_dot * omega**2
