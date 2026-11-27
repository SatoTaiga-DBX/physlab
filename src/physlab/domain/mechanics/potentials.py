"""衝突・散乱の実験で使う中心力ポテンシャル ``U(r)``（無限遠で 0）。

設計書: docs/design/007-collisions.md ／ 教科書: 第 5 章

どれも ``value(r)`` と ``derivative(r) = dU/dr`` を持つ。``breakpoint_radius`` は
式が切り替わる半径（求積の区間をそこで分ける）。なければ ``None``。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


class CentralPotential(Protocol):
    breakpoint_radius: float | None

    def value(self, r: Array) -> Array: ...

    def derivative(self, r: Array) -> Array: ...


@dataclass(frozen=True)
class PowerLaw:
    """``U = α / rⁿ``。``α > 0`` が斥力。``n = 1`` はクーロン力・万有引力の形。"""

    alpha: float = 1.0
    n: float = 1.0
    breakpoint_radius: float | None = None

    def value(self, r: Array) -> Array:
        return self.alpha / r**self.n

    def derivative(self, r: Array) -> Array:
        return -self.n * self.alpha / r ** (self.n + 1.0)


def coulomb(alpha: float = 1.0) -> PowerLaw:
    """点電荷どうしの ``U = α / r``（``α = k q Q``）。"""
    return PowerLaw(alpha=alpha, n=1.0)


@dataclass(frozen=True)
class Gaussian:
    """なめらかな丘（``height < 0`` なら谷） ``U = U0 exp(−r²/a²)``。"""

    height: float = 1.0
    width: float = 1.0
    breakpoint_radius: float | None = None

    def value(self, r: Array) -> Array:
        return self.height * np.exp(-((r / self.width) ** 2))

    def derivative(self, r: Array) -> Array:
        return -2.0 * r / self.width**2 * self.value(r)


@dataclass(frozen=True)
class UniformBall:
    """半径 R の球に一様に広がった電荷がつくる斥力。

    外（``r ≥ R``）は点電荷と同じ ``α/r``。中は、内側にある電荷だけが力を及ぼすので
    ``α (3 − r²/R²) / (2R)``。球の表面で値と傾きがつながる。
    """

    alpha: float = 1.0
    radius: float = 1.0

    @property
    def breakpoint_radius(self) -> float:
        return self.radius

    def value(self, r: Array) -> Array:
        inside = self.alpha * (3.0 - (r / self.radius) ** 2) / (2.0 * self.radius)
        outside = self.alpha / np.maximum(r, self.radius)
        return np.where(r < self.radius, inside, outside)

    def derivative(self, r: Array) -> Array:
        inside = -self.alpha * r / self.radius**3
        outside = -self.alpha / np.maximum(r, self.radius) ** 2
        return np.where(r < self.radius, inside, outside)


@dataclass(frozen=True)
class Yukawa:
    """遮蔽された斥力 ``α e^{−r/λ} / r``。"""

    alpha: float = 1.0
    length: float = 1.0
    breakpoint_radius: float | None = None

    def value(self, r: Array) -> Array:
        return self.alpha * np.exp(-r / self.length) / r

    def derivative(self, r: Array) -> Array:
        return -self.value(r) * (1.0 / r + 1.0 / self.length)


@dataclass(frozen=True)
class SoftDisk:
    """重なった分に比例して押し返すばね ``½ k (a − r)²``（``r < a``）、外では 0。"""

    stiffness: float = 4.0
    radius: float = 1.0

    @property
    def breakpoint_radius(self) -> float:
        return self.radius

    def value(self, r: Array) -> Array:
        return 0.5 * self.stiffness * np.maximum(self.radius - r, 0.0) ** 2

    def derivative(self, r: Array) -> Array:
        return -self.stiffness * np.maximum(self.radius - r, 0.0)
