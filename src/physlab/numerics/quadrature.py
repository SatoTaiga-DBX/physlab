"""ガウス・ルジャンドル求積。端点に 1/√ の特異性がある積分も扱う（学習用に自前で小さく実装）。

積分区間の端点 ``a``, ``b`` は配列でもよい（``(M,)``）。そのときは被積分関数に ``(M, n)`` の点が
渡るので、パラメータは ``x[:, None]`` のようにして一緒に並べる。
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray

Array = NDArray[np.float64]


def gauss_legendre(n: int) -> tuple[Array, Array]:
    """区間 [0, 1] の ``n`` 点ガウス・ルジャンドル求積の節点と重み。"""
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


def integrate_smooth(
    f: Callable[[Array], Array], a: ArrayLike, b: ArrayLike, n: int = 120
) -> Array:
    """``∫_a^b f(x) dx``。``f`` はなめらかで、端点でも有限とする。"""
    nodes, weights = gauss_legendre(n)
    a = np.asarray(a, dtype=float)[..., None]
    b = np.asarray(b, dtype=float)[..., None]
    x = a + (b - a) * nodes
    return np.sum(f(x) * weights, axis=-1) * (b - a)[..., 0]


def integrate_inverse_sqrt_end(
    f: Callable[[Array], Array], a: ArrayLike, b: ArrayLike, n: int = 120
) -> Array:
    """``∫_a^b f(x) dx``。``x → b`` で ``f ~ 1/√(b − x)`` の特異性があってもよい。

    ``x = b − (b − a) s²`` と置換すると ``dx = 2 (b − a) s ds`` の ``s`` が特異性を打ち消し、
    被積分関数がなめらかになる。端点そのものは評価しない。
    """
    nodes, weights = gauss_legendre(n)
    a = np.asarray(a, dtype=float)[..., None]
    b = np.asarray(b, dtype=float)[..., None]
    x = b - (b - a) * nodes**2
    return np.sum(f(x) * 2.0 * nodes * weights, axis=-1) * (b - a)[..., 0]
