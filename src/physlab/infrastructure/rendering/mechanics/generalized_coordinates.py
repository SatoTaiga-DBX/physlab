"""一般座標の図と動画（同じ位置・違う速度、二重振り子の配位空間）。

設計書: docs/design/005-analytical-mechanics-ch1.md（§3.1, §5）
結果データ（application）を受け取って描くだけで、計算はしない。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from physlab.application.mechanics.generalized_coordinates import (
    DoublePendulumRun,
    SameCoordinateRuns,
)
from physlab.domain.mechanics.double_pendulum import positions
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_animation, save_figure


def render_same_coordinate_figure(runs: SameCoordinateRuns, path: Path, t_max: float = 4.0) -> Path:
    """同じ ``θ0`` から出発しても、初期角速度が違えば運動は違う（状態は ``(q, q̇)`` で決まる）。"""
    setup_matplotlib()
    fig, (ax_t, ax_p) = plt.subplots(1, 2, figsize=(6.4, 2.9), constrained_layout=True)
    theta0_deg = np.rad2deg(runs.theta0)
    shown = runs.t <= t_max
    for k, w in enumerate(runs.omega0s):
        color = f"C{k}"
        label = f"$\\dot\\theta(0)={w:+.0f}$ rad/s"
        ax_t.plot(runs.t[shown], np.rad2deg(runs.theta[k][shown]), color=color, lw=1.8, label=label)
        ax_p.plot(np.rad2deg(runs.theta[k]), runs.omega[k], color=color, lw=1.4)
        ax_p.plot([theta0_deg], [w], "o", color=color, ms=7, mec="k", zorder=5)
    ax_t.plot([0], [theta0_deg], "ko", ms=6, zorder=6)
    ax_t.set_xlabel("時間 $t$ [s]")
    ax_t.set_ylabel("角度 $\\theta$ [deg]")
    ax_t.set_title(f"同じ $\\theta(0)={theta0_deg:.0f}^\\circ$ から出発", fontsize=10)
    low, high = ax_t.get_ylim()
    ax_t.set_ylim(low, high + 0.45 * (high - low))  # 凡例を置く余白
    ax_t.legend(fontsize=7, loc="upper center", ncol=2)
    ax_p.axvline(theta0_deg, color="0.5", ls=":", lw=1)
    ax_p.set_xlabel("角度 $\\theta$ [deg]")
    ax_p.set_ylabel("角速度 $\\dot\\theta$ [rad/s]")
    ax_p.set_title("位相平面：出発点は縦一列", fontsize=10)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def _time_colored_line(ax, x: np.ndarray, y: np.ndarray, t: np.ndarray, **kwargs) -> LineCollection:
    points = np.column_stack([x, y]).reshape(-1, 1, 2)
    collection = LineCollection(np.concatenate([points[:-1], points[1:]], axis=1), **kwargs)
    collection.set_array(t[:-1])
    ax.add_collection(collection)
    return collection


def render_double_pendulum_figure(run: DoublePendulumRun, path: Path) -> Path:
    """左: 実空間の軌跡と姿勢。右: 配位空間 ``(θ1, θ2)`` での軌跡（1 点が系全体を表す）。"""
    setup_matplotlib()
    p = run.params
    x1, y1, x2, y2 = positions(run.theta1, run.theta2, p)
    fig, (ax_r, ax_c) = plt.subplots(1, 2, figsize=(6.4, 3.2), constrained_layout=True)

    # 実空間: 先端の軌跡（時刻で色づけ）と、3 つの時刻の姿勢
    lc = _time_colored_line(ax_r, x2, y2, run.t, cmap="viridis", lw=0.9)
    for frac, alpha in ((0.0, 1.0), (0.12, 0.55), (0.3, 0.3)):
        i = int(frac * (run.t.size - 1))
        ax_r.plot([0, x1[i], x2[i]], [0, y1[i], y2[i]], "-o", color="k", alpha=alpha, lw=1.6, ms=5)
    ax_r.plot([0], [0], "k^", ms=7)
    reach = p.l1 + p.l2
    ax_r.set_xlim(-1.05 * reach, 1.05 * reach)
    ax_r.set_ylim(-1.05 * reach, 1.05 * reach)
    ax_r.set_aspect("equal")
    ax_r.set_xlabel("$x$ [m]")
    ax_r.set_ylabel("$y$ [m]")
    ax_r.set_title("実空間：先端の軌跡", fontsize=10)

    # 配位空間: 点が 1 つ動く
    _time_colored_line(ax_c, run.theta1, run.theta2, run.t, cmap="viridis", lw=0.9)
    ax_c.plot([run.theta1[0]], [run.theta2[0]], "o", color="tab:red", ms=7, mec="k", zorder=5)
    pad = 0.1
    ax_c.set_xlim(run.theta1.min() - pad, run.theta1.max() + pad)
    ax_c.set_ylim(run.theta2.min() - pad, run.theta2.max() + pad)
    ax_c.set_xlabel("$\\theta_1$ [rad]")
    ax_c.set_ylabel("$\\theta_2$ [rad]")
    ax_c.set_title("配位空間：1 点が系全体を表す", fontsize=10)
    fig.colorbar(lc, ax=ax_c, shrink=0.8, label="時刻 $t$ [s]")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_double_pendulum_video(run: DoublePendulumRun, path: Path, fps: int = 30) -> Path:
    """二重振り子の動画。左は実空間、右は配位空間 ``(θ1, θ2)``。再生は実時間と一致する。"""
    setup_matplotlib()
    p = run.params
    x1, y1, x2, y2 = positions(run.theta1, run.theta2, p)
    dt = float(run.t[1] - run.t[0])
    frames = np.round(np.arange(0.0, run.t[-1], 1.0 / fps) / dt).astype(int)
    trail = int(2.5 / dt)  # 先端の軌跡は直近 2.5 秒だけ残す

    fig, (ax_r, ax_c) = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={"width_ratios": [1, 1]})
    reach = p.l1 + p.l2
    ax_r.set_xlim(-1.1 * reach, 1.1 * reach)
    ax_r.set_ylim(-1.1 * reach, 1.1 * reach)
    ax_r.set_aspect("equal")
    ax_r.set_title("実空間（棒でつながった 2 つの錘）")
    ax_r.grid(False)
    ax_r.plot([0], [0], "k^", ms=9)
    (tail,) = ax_r.plot([], [], "-", color="tab:orange", lw=1.2, alpha=0.8)
    (rods,) = ax_r.plot([], [], "-", color="k", lw=2.5)
    (bob1,) = ax_r.plot([], [], "o", color="tab:blue", ms=16)
    (bob2,) = ax_r.plot([], [], "o", color="tab:red", ms=16)
    time_text = ax_r.text(0.03, 0.04, "", transform=ax_r.transAxes)

    ax_c.set_xlim(run.theta1.min() - 0.15, run.theta1.max() + 0.15)
    ax_c.set_ylim(run.theta2.min() - 0.15, run.theta2.max() + 0.15)
    ax_c.set_xlabel("$\\theta_1$ [rad]")
    ax_c.set_ylabel("$\\theta_2$ [rad]")
    ax_c.set_title("配位空間 $(\\theta_1,\\theta_2)$：1 点が系全体の状態")
    ax_c.plot(run.theta1, run.theta2, color="0.85", lw=0.8)
    (trace,) = ax_c.plot([], [], color="tab:purple", lw=1.2)
    (point,) = ax_c.plot([], [], "o", color="tab:purple", ms=10)
    fig.tight_layout()

    def update(i: int):
        lo = max(0, i - trail)
        tail.set_data(x2[lo : i + 1], y2[lo : i + 1])
        rods.set_data([0, x1[i], x2[i]], [0, y1[i], y2[i]])
        bob1.set_data([x1[i]], [y1[i]])
        bob2.set_data([x2[i]], [y2[i]])
        trace.set_data(run.theta1[: i + 1], run.theta2[: i + 1])
        point.set_data([run.theta1[i]], [run.theta2[i]])
        time_text.set_text(f"t = {run.t[i]:5.2f} s")
        return tail, rods, bob1, bob2, trace, point, time_text

    anim = animation.FuncAnimation(fig, update, frames=frames, blit=True)
    saved = save_animation(anim, path, fps=fps)
    plt.close(fig)
    return saved
