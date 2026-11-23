"""1 次元の質点のラグランジアン・作用・停留経路（最小作用の原理の道具立て）。

設計書: docs/design/005-analytical-mechanics-ch1.md ／ 教科書: 第 3 章「運動方程式」§最小作用の原理

ラグランジアン ``L(q, q̇) = ½ m q̇² − V(q)`` と作用 ``S[q] = ∫ L dt``。
両端 ``(t1, q1)``, ``(t2, q2)`` を固定したとき、真の運動は ``S`` の停留点になる。
ここでは 3 つの系（自由粒子・一様重力・調和振動子）について、停留経路を解析的に持つ。

3 つの系はいずれも ``L`` が ``(q, q̇)`` の 2 次式なので、停留経路 ``q*`` に
``ε η_n``（``η_n = sin(nπ(t-t1)/τ)``、両端で 0）を足した経路の作用は ``ε`` のちょうど 2 次式になる::

    S(ε) = S* + c_n ε²,    c_n = (m τ / 4) [ (nπ/τ)² − ω² ]

最内層（domain）なので、純粋な関数とデータだけを置く。
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import simpson

Array = NDArray[np.float64]
# 停留経路を返す関数 path(t, t1, q1, t2, q2)
StationaryPath = Callable[[Array, float, float, float, float], Array]


@dataclass(frozen=True)
class System1D:
    """``L = ½ m q̇² − V(q)`` で表される 1 次元の質点。"""

    name: str
    mass: float
    omega: float  # 調和振動子の角振動数（自由・一様重力は 0）。2 次の係数の公式に使う
    potential: Callable[[Array], Array]  # V(q)
    dpotential: Callable[[Array], Array]  # V'(q)
    stationary: StationaryPath  # オイラー＝ラグランジュ方程式の解（端点つき）

    def lagrangian(self, q: Array, qdot: Array) -> Array:
        """``L = ½ m q̇² − V(q)``。"""
        return 0.5 * self.mass * qdot**2 - self.potential(q)

    def stationary_path(self, t: Array, t1: float, q1: float, t2: float, q2: float) -> Array:
        """両端を ``(t1, q1)``, ``(t2, q2)`` に固定したときの真の運動 ``q*(t)``。"""
        return self.stationary(t, t1, q1, t2, q2)

    def action_curvature(self, tau: float, n: int = 1) -> float:
        """``S(ε) = S* + c_n ε²`` の ``c_n``。正なら ``η_n`` の方向では最小、負なら極大。"""
        return 0.25 * self.mass * tau * ((n * np.pi / tau) ** 2 - self.omega**2)


def _free_path(t: Array, t1: float, q1: float, t2: float, q2: float) -> Array:
    return q1 + (q2 - q1) * (t - t1) / (t2 - t1)


def free_particle(mass: float = 1.0) -> System1D:
    """自由粒子 ``V = 0``。真の運動は両端を結ぶ直線（等速）。"""
    return System1D(
        "自由粒子", mass, 0.0, lambda q: np.zeros_like(q), lambda q: np.zeros_like(q), _free_path
    )


def uniform_gravity(mass: float = 1.0, g: float = 9.81) -> System1D:
    """一様重力 ``V = m g q``（``q`` は鉛直上向きの高さ）。真の運動は放物線。"""

    def path(t: Array, t1: float, q1: float, t2: float, q2: float) -> Array:
        s, tau = t - t1, t2 - t1
        return q1 + (q2 - q1) * s / tau + 0.5 * g * s * (tau - s)

    return System1D(
        "一様重力", mass, 0.0, lambda q: mass * g * q, lambda q: np.full_like(q, mass * g), path
    )


def harmonic_oscillator(mass: float = 1.0, omega: float = 2.0 * np.pi) -> System1D:
    """調和振動子 ``V = ½ m ω² q²``。真の運動は ``sin ωτ ≠ 0`` のとき一意に決まる。"""

    def path(t: Array, t1: float, q1: float, t2: float, q2: float) -> Array:
        s, tau = t - t1, t2 - t1
        denominator = np.sin(omega * tau)
        if abs(denominator) < 1e-9:
            raise ValueError("sin(ωτ) = 0: 両端を結ぶ真の運動が一意に決まらない（τ = nπ/ω）")
        return (q1 * np.sin(omega * (tau - s)) + q2 * np.sin(omega * s)) / denominator

    return System1D(
        "調和振動子",
        mass,
        omega,
        lambda q: 0.5 * mass * omega**2 * q**2,
        lambda q: mass * omega**2 * q,
        path,
    )


def sine_mode(t: Array, t1: float, t2: float, n: int = 1) -> Array:
    """両端で 0 になる変分の形 ``η_n(t) = sin(nπ(t-t1)/τ)``。"""
    return np.sin(n * np.pi * (t - t1) / (t2 - t1))


# ---- 連続の作用 ------------------------------------------------------------------------


def action(system: System1D, t: Array, q: Array, qdot: Array | None = None) -> float:
    """作用 ``S = ∫ L dt``（シンプソン則）。``qdot`` を省略すると ``q`` の数値微分を使う。"""
    if qdot is None:
        qdot = np.gradient(q, t, edge_order=2)
    return float(simpson(system.lagrangian(q, qdot), x=t))


def action_batch(system: System1D, t: Array, paths: Array) -> Array:
    """複数の経路（形状 ``(n, len(t))``）の作用を一度に求める。"""
    qdot = np.gradient(paths, t, axis=-1, edge_order=2)
    return np.asarray(simpson(system.lagrangian(paths, qdot), x=t, axis=-1))


# ---- 離散の作用 ------------------------------------------------------------------------


def discrete_action(system: System1D, t: Array, q: Array) -> float:
    """離散作用 ``Σ Δt [½ m ((q_{k+1}-q_k)/Δt)² − ½ (V_k + V_{k+1})]``（``t`` は等間隔）。"""
    dt = t[1] - t[0]
    velocity = np.diff(q) / dt
    v = system.potential(q)
    return float(np.sum(dt * (0.5 * system.mass * velocity**2 - 0.5 * (v[:-1] + v[1:]))))


def discrete_action_gradient(system: System1D, t: Array, q: Array) -> Array:
    """離散作用の各節点での勾配。両端は固定なので 0 を返す。

    内部の節点では ``g_k = (m/Δt)(2 q_k − q_{k-1} − q_{k+1}) − Δt V'(q_k)``。
    ``g_k = 0`` は ``m (q_{k+1} − 2 q_k + q_{k-1}) / Δt² = −V'(q_k)``（速度ベルレ法の位置の式）。
    """
    dt = t[1] - t[0]
    grad = np.zeros_like(q)
    grad[1:-1] = system.mass / dt * (2.0 * q[1:-1] - q[:-2] - q[2:]) - dt * system.dpotential(
        q[1:-1]
    )
    return grad
