"""崩壊の実験: 乱数で向きを選んで何度も崩壊させ、実験室系の角度とエネルギーを数える。

設計書: docs/design/007-collisions.md ／ 教科書: 第 5 章 §粒子の崩壊

実験室系の量は、3 次元の速度ベクトルを足し合わせて直接求める（角度分布の公式は使わない）。
理論の式（``domain/mechanics/decay.py``）との比較は、数えた結果に対して行う。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from physlab.domain.mechanics.decay import (
    DecayParams,
    angle_density,
    kinetic_energy_range,
    max_lab_angle,
)
from physlab.numerics.quadrature import Array, integrate_inverse_sqrt_end


@dataclass(frozen=True)
class DecaySample:
    """多数回の崩壊の結果。親は x 軸の向きに速さ ``v_parent`` で動いている。"""

    params: DecayParams
    v_parent: float
    v1: Array  # (n, 3) 破片 1 の実験室系の速度
    v2: Array  # (n, 3) 破片 2
    theta: Array  # (n,) 破片 1 の向き（親の進行方向から）
    kinetic: Array  # (n,) 破片 1 の運動エネルギー

    @property
    def speed_ratio(self) -> float:
        """``k = V / v10``。破片 1 の角度分布を決める。"""
        return self.v_parent / self.params.rest_speeds[0]


@dataclass(frozen=True)
class Histogram:
    """区間ごとの確率（実測と理論）。"""

    edges: Array
    measured: Array
    expected: Array
    error: Array  # 実測の統計誤差（1σ）

    @property
    def centers(self) -> Array:
        return 0.5 * (self.edges[:-1] + self.edges[1:])

    @property
    def widths(self) -> Array:
        return np.diff(self.edges)


def isotropic_directions(n: int, rng: np.random.Generator) -> Array:
    """球面上に一様な単位ベクトル ``(n, 3)``。x 軸からの角の余弦が一様分布。"""
    cos_t = rng.uniform(-1.0, 1.0, n)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    sin_t = np.sqrt(1.0 - cos_t**2)
    return np.column_stack([cos_t, sin_t * np.cos(phi), sin_t * np.sin(phi)])


def simulate_decays(
    params: DecayParams, v_parent: float, n: int = 200_000, seed: int = 1
) -> DecaySample:
    direction = isotropic_directions(n, np.random.default_rng(seed))
    drift = np.array([v_parent, 0.0, 0.0])
    p0 = params.momentum
    v1 = drift + p0 / params.m1 * direction
    v2 = drift - p0 / params.m2 * direction
    speed = np.linalg.norm(v1, axis=1)
    theta = np.arccos(np.clip(v1[:, 0] / speed, -1.0, 1.0))
    return DecaySample(params, v_parent, v1, v2, theta, 0.5 * params.m1 * speed**2)


def momentum_error(sample: DecaySample) -> float:
    """``m1 v1 + m2 v2 − (m1 + m2) V`` の最大の大きさ。"""
    p = sample.params
    total = p.m1 * sample.v1 + p.m2 * sample.v2
    total[:, 0] -= (p.m1 + p.m2) * sample.v_parent
    return float(np.max(np.linalg.norm(total, axis=1)))


def energy_error(sample: DecaySample) -> float:
    """破片の運動エネルギーの和と、``½ M V² + ε`` の差の最大値。"""
    p = sample.params
    after = 0.5 * p.m1 * np.sum(sample.v1**2, axis=1) + 0.5 * p.m2 * np.sum(sample.v2**2, axis=1)
    before = 0.5 * (p.m1 + p.m2) * sample.v_parent**2 + p.energy
    return float(np.max(np.abs(after - before)))


def angle_bin_probability(edges: Array, k: float) -> Array:
    """理論の密度を区間ごとに積分した確率。上限の角の特異性は置換で打ち消す。"""
    limit = max_lab_angle(k)
    low = np.minimum(edges[:-1], limit)
    high = np.minimum(edges[1:], limit)
    probability = np.zeros(low.size)
    inside = high > low
    probability[inside] = integrate_inverse_sqrt_end(
        lambda t: angle_density(t, k), low[inside], high[inside], n=200
    )
    return probability


def _histogram(values: Array, edges: Array, expected: Array) -> Histogram:
    counts, _ = np.histogram(values, bins=edges)
    n = values.size
    return Histogram(edges, counts / n, expected, np.sqrt(np.maximum(counts, 1)) / n)


def angle_histogram(sample: DecaySample, n_bins: int = 36) -> Histogram:
    """実験室系の角の分布（0 から上限の角まで）。"""
    k = sample.speed_ratio
    edges = np.linspace(0.0, max_lab_angle(k), n_bins + 1)
    edges[-1] += 1e-9  # 上限ちょうどの値を最後の区間に含める
    return _histogram(sample.theta, edges, angle_bin_probability(edges, k))


def energy_histogram(sample: DecaySample, n_bins: int = 20) -> Histogram:
    """破片 1 の運動エネルギーの分布。理論は一様。"""
    low, high = kinetic_energy_range(
        sample.params.m1, sample.v_parent, sample.params.rest_speeds[0]
    )
    edges = np.linspace(low, high, n_bins + 1)
    edges[0] -= 1e-9 * (high - low)
    edges[-1] += 1e-9 * (high - low)
    return _histogram(sample.kinetic, edges, np.full(n_bins, 1.0 / n_bins))


def planar_burst(params: DecayParams, n: int = 72) -> Array:
    """親の静止系で、1 つの平面内に等間隔の向きへ飛ぶ破片 1 の速度 ``(n, 2)``（動画用）。"""
    angle = 2.0 * np.pi * np.arange(n) / n
    return params.rest_speeds[0] * np.column_stack([np.cos(angle), np.sin(angle)])


@dataclass(frozen=True)
class Burst:
    """x 軸を進む親が、時刻 0・原点で崩壊する様子。"""

    v_parent: float
    v0: float
    t: Array
    parent_x: Array  # (len(t),) 親（崩壊後は静止系の原点）の位置
    fragments: Array  # (len(t), n, 2)。崩壊前（t < 0）は親の位置に重なる


def decay_burst(params: DecayParams, v_parent: float, t: Array, n: int = 72) -> Burst:
    velocity = planar_burst(params, n) + np.array([v_parent, 0.0])
    parent = v_parent * t
    after = t[:, None, None] * velocity[None, :, :]
    before = np.zeros_like(after)
    before[:, :, 0] = parent[:, None]
    fragments = np.where((t < 0)[:, None, None], before, after)
    return Burst(v_parent, params.rest_speeds[0], t, parent, fragments)
