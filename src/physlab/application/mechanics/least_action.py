"""最小作用の原理のユースケース: 経路の族の作用、作用の曲率、作用の地形、経路の緩和。

設計書: docs/design/005-analytical-mechanics-ch1.md ／ 教科書: 第 3 章「運動方程式」§最小作用の原理
返すのは dataclass（配列と派生量）だけ。これを描画するのは infrastructure の仕事。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from physlab.domain.mechanics.lagrangian import (
    System1D,
    action,
    action_batch,
    discrete_action,
    discrete_action_gradient,
    sine_mode,
)
from physlab.numerics.integrators import Array
from physlab.numerics.optimize import steepest_descent

# ---- 経路の族と作用 ----------------------------------------------------------------------


@dataclass(frozen=True)
class ActionScan:
    """真の運動 ``q*`` に ``ε η_n`` を足した経路の族と、その作用 ``S(ε)``。"""

    system: System1D
    t: Array
    true_path: Array
    mode: int
    eps: Array
    paths: Array  # 形状 (len(eps), len(t))
    action: Array  # S(ε)
    action_true: float  # S*（真の運動の作用）
    fit_curvature: float  # S(ε) を 2 次式で当てはめた ε² の係数 c
    fit_linear: float  # 同 ε の係数（停留なら 0）
    theory_curvature: float  # 公式 (mτ/4)[(nπ/τ)² − ω²]


def scan_action(
    system: System1D,
    t1: float,
    q1: float,
    t2: float,
    q2: float,
    eps: Array,
    mode: int = 1,
    n_points: int = 4001,
) -> ActionScan:
    """両端 ``(t1, q1)``, ``(t2, q2)`` を固定した経路の族の作用を求める。"""
    t = np.linspace(t1, t2, n_points)
    q_true = system.stationary_path(t, t1, q1, t2, q2)
    eps = np.asarray(eps, dtype=float)
    paths = q_true[None, :] + eps[:, None] * sine_mode(t, t1, t2, mode)[None, :]
    s = action_batch(system, t, paths)
    c, b, _ = np.polyfit(eps, s, 2)
    return ActionScan(
        system=system,
        t=t,
        true_path=q_true,
        mode=mode,
        eps=eps,
        paths=paths,
        action=s,
        action_true=action(system, t, q_true),
        fit_curvature=float(c),
        fit_linear=float(b),
        theory_curvature=system.action_curvature(t2 - t1, mode),
    )


# ---- 作用の曲率（最小か鞍点か） -----------------------------------------------------------


@dataclass(frozen=True)
class CurvatureScan:
    """時間の長さ ``τ`` を変えたときの、各モードの作用の曲率 ``c_n(τ)``。"""

    system: System1D
    tau: Array
    modes: tuple[int, ...]
    numeric: dict[int, Array]  # 数値積分と中心差分から求めた c_n
    theory: dict[int, Array]  # 公式


def scan_curvature_vs_duration(
    system: System1D,
    taus: Array,
    modes: tuple[int, ...] = (1, 2, 3),
    eps: float = 0.1,
    n_points: int = 2001,
) -> CurvatureScan:
    """``c_n(τ)`` を、ゼロの経路のまわりの ``ε η_n`` で数値的に求める。

    ``S`` は ``ε`` の 2 次式なので、中心差分 ``(S(ε)+S(−ε)−2S(0)) / (2ε²)`` が ``c_n`` に一致する。
    ゼロの経路（両端 0）は ``V(0) = 0`` の系で ``τ = nπ/ω`` でも定義できる。
    """
    taus = np.asarray(taus, dtype=float)
    numeric = {n: np.empty_like(taus) for n in modes}
    theory = {n: np.empty_like(taus) for n in modes}
    for i, tau in enumerate(taus):
        t = np.linspace(0.0, tau, n_points)
        for n in modes:
            eta = sine_mode(t, 0.0, tau, n)
            s_plus = action(system, t, eps * eta)
            s_minus = action(system, t, -eps * eta)
            s_zero = action(system, t, np.zeros_like(t))
            numeric[n][i] = (s_plus + s_minus - 2.0 * s_zero) / (2.0 * eps**2)
            theory[n][i] = system.action_curvature(tau, n)
    return CurvatureScan(system, taus, tuple(modes), numeric, theory)


# ---- 作用の地形 --------------------------------------------------------------------------


@dataclass(frozen=True)
class ActionLandscape:
    """2 つのモード ``a η_1 + b η_2`` の張る平面での ``S(a, b) − S*``。"""

    system: System1D
    tau: float
    a: Array
    b: Array
    delta_action: Array  # 形状 (len(b), len(a))
    curvature_a: float  # 理論値 c_1
    curvature_b: float  # 理論値 c_2


def action_landscape(
    system: System1D,
    t1: float,
    q1: float,
    t2: float,
    q2: float,
    half_width: float = 1.0,
    n: int = 41,
    n_points: int = 801,
) -> ActionLandscape:
    t = np.linspace(t1, t2, n_points)
    q_true = system.stationary_path(t, t1, q1, t2, q2)
    a = np.linspace(-half_width, half_width, n)
    b = np.linspace(-half_width, half_width, n)
    grid_a, grid_b = np.meshgrid(a, b)
    paths = (
        q_true[None, :]
        + grid_a.reshape(-1, 1) * sine_mode(t, t1, t2, 1)[None, :]
        + grid_b.reshape(-1, 1) * sine_mode(t, t1, t2, 2)[None, :]
    )
    s = action_batch(system, t, paths).reshape(grid_a.shape)
    tau = t2 - t1
    return ActionLandscape(
        system=system,
        tau=tau,
        a=a,
        b=b,
        delta_action=s - action(system, t, q_true),
        curvature_a=system.action_curvature(tau, 1),
        curvature_b=system.action_curvature(tau, 2),
    )


# ---- 経路の緩和（作用を下げ続けると真の運動に収束する） ---------------------------------------


@dataclass(frozen=True)
class PathRelaxation:
    """離散作用を最急降下法で下げていく途中の経路。"""

    system: System1D
    t: Array  # 節点の時刻
    true_path: Array  # 解析的な停留経路（節点上）
    iterations: Array  # 記録した反復回数（先頭は 0）
    paths: Array  # 形状 (len(iterations), len(t))
    action: Array  # 各時点の離散作用
    max_error: Array  # max |q_k − q*_k|
    action_true: float  # 解析的な停留経路の離散作用


def wiggly_initial_path(t: Array, t1: float, q1: float, t2: float, q2: float) -> Array:
    """両端を結ぶ直線に、粗い波と細かい波を足した、決定的な初期経路。"""
    s = (t - t1) / (t2 - t1)
    return q1 + (q2 - q1) * s + 2.5 * np.sin(3 * np.pi * s) + 0.9 * np.sin(11 * np.pi * s)


def relax_path(
    system: System1D,
    t1: float,
    q1: float,
    t2: float,
    q2: float,
    initial: Array | None = None,
    n_nodes: int = 41,
    n_iter: int = 8000,
    n_frames: int = 96,
    frame_horizon: int = 3000,
    safety: float = 0.9,
) -> PathRelaxation:
    """離散作用を最急降下法で下げ、途中の経路を ``n_frames`` 個前後（対数間隔）記録する。

    刻み ``α = safety · Δt / (2m)`` は、勾配の最大固有値 ``4m/Δt`` から決まる
    安定性の上限の ``safety`` 倍。
    停留点が鞍点（``c_1 ≤ 0``）の系では最小化できないので ``ValueError``。
    """
    tau = t2 - t1
    if system.action_curvature(tau, 1) <= 0.0:
        raise ValueError(
            f"τ = {tau:g} では停留点が鞍点（c_1 ≤ 0）で、作用を下げても真の運動に収束しない"
        )
    t = np.linspace(t1, t2, n_nodes)
    q0 = wiggly_initial_path(t, t1, q1, t2, q2) if initial is None else np.array(initial, float)
    q0[0], q0[-1] = q1, q2  # 両端は固定
    step = safety * (t[1] - t[0]) / (2.0 * system.mass)

    horizon = min(frame_horizon, n_iter)
    record = {0, n_iter, *np.round(np.geomspace(1, horizon, n_frames - 2)).astype(int).tolist()}
    _, snapshots = steepest_descent(
        lambda q: discrete_action_gradient(system, t, q), q0, step, n_iter, record
    )
    iterations = np.array(sorted(snapshots))
    paths = np.array([snapshots[k] for k in iterations])
    true_path = system.stationary_path(t, t1, q1, t2, q2)
    return PathRelaxation(
        system=system,
        t=t,
        true_path=true_path,
        iterations=iterations,
        paths=paths,
        action=np.array([discrete_action(system, t, p) for p in paths]),
        max_error=np.max(np.abs(paths - true_path[None, :]), axis=1),
        action_true=discrete_action(system, t, true_path),
    )
