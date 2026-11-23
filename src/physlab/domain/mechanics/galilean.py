"""ガリレイ変換と、自由粒子の作用のずれ。

設計書: docs/design/005-analytical-mechanics-ch1.md ／ 教科書: 第 3 章 §ガリレイの相対性原理

慣性系 K' が K に対して速度 V で動くとき（ランダウ＝リフシッツの向き）::

    r = r' + V t,   v = v' + V,   t = t'

自由粒子のラグランジアンは ``L = ½ m v² = L' + df/dt``、``f = m V·r' + ½ m V² t``。
よって両端を固定した作用は ``S − S' = f(t2) − f(t1)`` と、**経路によらない定数**だけずれる。
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


def to_moving_frame_position(r: Array, t: Array, v_frame: Array | float) -> Array:
    """K での位置 ``r(t)`` を、速度 ``V`` で動く K' での位置 ``r' = r − V t`` にする。

    ``r`` は形状 ``(n,)``（1 次元）か ``(n, d)``、``V`` はスカラーか長さ ``d`` のベクトル。
    """
    t = np.asarray(t, dtype=float)
    v_frame = np.asarray(v_frame, dtype=float)
    return r - (t * v_frame if r.ndim == 1 else t[:, None] * v_frame)


def to_moving_frame_velocity(v: Array, v_frame: Array | float) -> Array:
    """K での速度を K' での速度 ``v' = v − V`` にする。"""
    return v - np.asarray(v_frame, dtype=float)


def free_particle_action_shift(mass: float, v_frame: float, dq_prime: float, tau: float) -> float:
    """1 次元の自由粒子で ``S − S' = m V Δq' + ½ m V² τ``（``Δq' = q'(t2) − q'(t1)``）。"""
    return mass * v_frame * dq_prime + 0.5 * mass * v_frame**2 * tau
