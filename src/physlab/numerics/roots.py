"""2 分法。配列でまとめて解く（物理を知らない）。"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray

Array = NDArray[np.float64]


def bisect(f: Callable[[Array], Array], lo: ArrayLike, hi: ArrayLike, n_iter: int = 100) -> Array:
    """``f(lo) > 0`` と ``f(hi) <= 0`` の間の零点を挟み込み、**正の側の端**を返す。

    正の側を返すので、``f`` の平方根を ``[lo, 返り値]`` で安全に取れる。
    ``lo``・``hi`` は同じ形の配列。``f(lo)`` は呼ばない（``lo`` で定義できない関数でもよい）。

    Raises:
        ValueError: ``f(hi) > 0`` の要素がある（零点が挟まれていない）。
    """
    low = np.array(lo, dtype=float)
    high = np.array(hi, dtype=float)
    if np.any(f(high) > 0):
        raise ValueError("f(hi) <= 0 でない区間がある（零点が挟まれていない）")
    for _ in range(n_iter):
        mid = 0.5 * (low + high)
        positive = f(mid) > 0
        low = np.where(positive, mid, low)
        high = np.where(positive, high, mid)
    return low
