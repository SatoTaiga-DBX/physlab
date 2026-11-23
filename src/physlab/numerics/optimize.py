"""物理を知らない最適化。最急降下法だけを置く（学習用に自前で小さく実装）。"""

from __future__ import annotations

from collections.abc import Callable, Iterable

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


def steepest_descent(
    grad: Callable[[Array], Array],
    x0: Array,
    step: float,
    n_iter: int,
    record_at: Iterable[int] = (),
) -> tuple[Array, dict[int, Array]]:
    """``x ← x − step · grad(x)`` を ``n_iter`` 回繰り返す。

    Args:
        record_at: その反復回数（0 は初期値）での ``x`` を途中経過として残す。

    Returns:
        ``(最終の x, {反復回数: その時点の x のコピー})``。
    """
    wanted = set(record_at)
    x = np.array(x0, dtype=float)
    snapshots: dict[int, Array] = {}
    if 0 in wanted:
        snapshots[0] = x.copy()
    for k in range(1, n_iter + 1):
        x = x - step * grad(x)
        if k in wanted:
            snapshots[k] = x.copy()
    return x, snapshots
