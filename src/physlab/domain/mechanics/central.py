"""中心力 ``U = −α / rⁿ`` の中の質点（2 次元）。角運動量と面積速度。

設計書: docs/design/006-conservation-laws.md ／ 教科書: 第 4 章 §角運動量

``0 < n < 2`` なら遠心力の障壁があるので、中心に落ち込まずに回り続ける（``n = 1`` は万有引力）。
``n ≠ 1`` の軌道は閉じず、ロゼット状に回転する。それでも面積速度 ``L / (2m)`` は一定。
状態: 位置 ``q = [x, y]``、速度 ``v = [vx, vy]``。
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


@dataclass(frozen=True)
class CentralParams:
    alpha: float = 1.0
    exponent: float = 1.4  # n（0 < n < 2）
    mass: float = 1.0

    def __post_init__(self) -> None:
        if not 0.0 < self.exponent < 2.0:
            raise ValueError("exponent は 0 < n < 2（中心に落ち込まない範囲）")


def accel(p: CentralParams) -> Callable[[float, Array], Array]:
    """``a = −α n r^{−n−2} q / m``。"""

    def a(_t: float, q: Array) -> Array:
        r = np.linalg.norm(q, axis=-1, keepdims=True)
        return -p.alpha * p.exponent * r ** (-p.exponent - 2.0) * q / p.mass

    return a


def energy(q: Array, v: Array, p: CentralParams) -> Array:
    r = np.linalg.norm(q, axis=-1)
    return 0.5 * p.mass * np.sum(v**2, axis=-1) - p.alpha / r**p.exponent


def angular_momentum(q: Array, v: Array, p: CentralParams) -> Array:
    return p.mass * (q[..., 0] * v[..., 1] - q[..., 1] * v[..., 0])


def circular_speed(radius: float, p: CentralParams) -> float:
    """半径 ``radius`` の円軌道の速さ ``√(α n / (m rⁿ))``。"""
    return float(np.sqrt(p.alpha * p.exponent / (p.mass * radius**p.exponent)))
