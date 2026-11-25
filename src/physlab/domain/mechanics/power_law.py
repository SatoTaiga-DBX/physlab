"""冪ポテンシャル ``U = |x|^k / k``（1 次元、m = 1）と力学的相似。

設計書: docs/design/006-conservation-laws.md ／ 教科書: 第 4 章 §力学的相似

振幅 A から静かに放すときの周期は、エネルギー保存から厳密に求まる::

    T(A) = 4 √(k/2) A^{1 − k/2} Γ(1/k) Γ(1/2) / (k Γ(1/k + 1/2))

相似則 ``x_A(t) = A x_1(t A^{k/2 − 1})`` から周期は ``A^{1 − k/2}`` に比例する。
ビリアル定理: 時間平均で ``2⟨T⟩ = k⟨U⟩``。
"""

from __future__ import annotations

from collections.abc import Callable
from math import gamma, pi, sqrt

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


def potential(x: Array, k: float) -> Array:
    return np.abs(x) ** k / k


def accel(k: float) -> Callable[[float, Array], Array]:
    """``a = −dU/dx = −sign(x) |x|^{k−1}``（``k ≥ 1``）。"""

    def a(_t: float, x: Array) -> Array:
        return -np.sign(x) * np.abs(x) ** (k - 1.0)

    return a


def period_exact(k: float, amplitude: float) -> float:
    """振幅 ``amplitude`` から静かに放したときの周期。"""
    integral = gamma(1.0 / k) * sqrt(pi) / (k * gamma(1.0 / k + 0.5))  # ∫0^1 du / √(1 − u^k)
    return float(4.0 * sqrt(k / 2.0) * amplitude ** (1.0 - k / 2.0) * integral)


def time_scale(k: float, amplitude: float) -> float:
    """相似則の時間の比。振幅 1 の運動を ``amplitude`` に写すと時間は ``A^{1 − k/2}`` 倍。"""
    return float(amplitude ** (1.0 - k / 2.0))
