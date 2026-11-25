"""面積速度の図と動画。結果データを受け取って描くだけで、計算はしない。

設計書: docs/design/006-conservation-laws.md
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon

from physlab.application.mechanics.areal_velocity import Sectors
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_animation, save_figure

SECTOR_COLORS = ("#4c72b0", "#dd8452")


def _wedge(sectors: Sectors, a: int, b: int) -> np.ndarray:
    """中心と、軌道の添字 a..b の点がつくる多角形。"""
    return np.vstack([[0.0, 0.0], sectors.orbit.q[a : b + 1]])


def _frame(ax, sectors: Sectors) -> None:
    q = sectors.orbit.q
    reach = 1.08 * np.abs(q).max()
    ax.set_xlim(-reach, reach)
    ax.set_ylim(-reach, reach)
    ax.set_aspect("equal")
    ax.plot([0], [0], "ko", ms=6)


def render_areal_figure(sectors: Sectors, path: Path) -> Path:
    """等しい時間の扇形。面積はすべて同じで、弧の長さは中心に近いほど長い。"""
    setup_matplotlib()
    fig = plt.figure(figsize=(6.4, 3.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[1.15, 1])
    ax = fig.add_subplot(grid[:, 0])
    ax_area = fig.add_subplot(grid[0, 1])
    ax_arc = fig.add_subplot(grid[1, 1])

    _frame(ax, sectors)
    q = sectors.orbit.q
    ax.plot(q[:, 0], q[:, 1], color="0.7", lw=0.8)
    edges = sectors.edges
    for k, (a, b) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        ax.add_patch(
            Polygon(
                _wedge(sectors, a, b), closed=True, alpha=0.55, color=SECTOR_COLORS[k % 2], lw=0
            )
        )
        ax.plot(*q[a], "k.", ms=4)
    ax.set_title(f"1 区間 = {sectors.duration:.3f} s の扇形", fontsize=10)
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")

    number = np.arange(1, sectors.areas.size + 1)
    colors = [SECTOR_COLORS[k % 2] for k in range(sectors.areas.size)]
    ax_area.bar(number, sectors.areas, color=colors)
    ax_area.set_ylim(0, 1.25 * sectors.areas.max())
    ax_area.set_title("掃いた面積 [m²]：すべて同じ", fontsize=9)
    ax_arc.bar(number, sectors.arc_lengths, color=colors)
    ax_arc.set_title("軌道に沿った長さ [m]：大きく違う", fontsize=9)
    ax_arc.set_xlabel("区間の番号")
    for axis in (ax_area, ax_arc):
        axis.set_xticks(number[::2])
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_areal_video(sectors: Sectors, path: Path, fps: int = 30) -> Path:
    """質点が動くにつれて扇形が塗られていく動画（再生は実時間）。"""
    setup_matplotlib()
    orbit = sectors.orbit
    q = orbit.q
    dt = float(orbit.t[1] - orbit.t[0])
    idx = np.round(np.arange(0.0, orbit.t[-1], 1.0 / fps) / dt).astype(int)
    steps = int(sectors.edges[1] - sectors.edges[0])

    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    _frame(ax, sectors)
    ax.set_title("同じ時間に掃く面積は同じ")
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")
    (trail,) = ax.plot([], [], color="0.35", lw=1)
    (radius,) = ax.plot([], [], color="k", lw=1.2)
    (body,) = ax.plot([], [], "o", color="tab:red", ms=11)
    text = ax.text(0.02, 0.96, "", transform=ax.transAxes, va="top")
    patches: list[Polygon] = []
    fig.tight_layout()

    def update(i: int):
        trail.set_data(q[: i + 1, 0], q[: i + 1, 1])
        radius.set_data([0, q[i, 0]], [0, q[i, 1]])
        body.set_data([q[i, 0]], [q[i, 1]])
        # 扇形は、区間の長さごとに色を変えて足していく（直近の 12 区間だけ残す）
        current = i // steps
        for patch in patches:
            patch.remove()
        patches.clear()
        for k in range(max(0, current - 11), current + 1):
            a, b = k * steps, min((k + 1) * steps, i)
            if b > a:
                patch = Polygon(
                    _wedge(sectors, a, b), closed=True, alpha=0.45, color=SECTOR_COLORS[k % 2], lw=0
                )
                ax.add_patch(patch)
                patches.append(patch)
        swept = orbit.step_area[:i].sum()
        text.set_text(f"t = {orbit.t[i]:.2f} s\n掃いた面積の合計 = {swept:.3f} m²")
        return (trail, radius, body, text, *patches)

    anim = animation.FuncAnimation(fig, update, frames=idx, blit=False)
    saved = save_animation(anim, path, fps=fps)
    plt.close(fig)
    return saved
