"""崩壊の図と動画。結果データを受け取って描くだけで、計算はしない。

設計書: docs/design/007-collisions.md
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt

from physlab.application.mechanics.decay import Burst, Histogram
from physlab.domain.mechanics.decay import angle_density, max_lab_angle
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_animation, save_figure

COLORS = ("tab:blue", "tab:red", "tab:green")


def _velocity_panel(ax, k: float) -> None:
    """速度の図: 親の速度 V の先を中心に、半径 v0 の円（v0 = 1 に規格化）。"""
    phi = np.linspace(0.0, 2.0 * np.pi, 400)
    ax.plot(k + np.cos(phi), np.sin(phi), color="0.4", lw=1.4)
    ax.annotate(
        "", xy=(k, 0), xytext=(0, 0), arrowprops={"arrowstyle": "->", "color": "k", "lw": 1.6}
    )
    ax.plot([0], [0], "ko", ms=4)
    ax.plot([k], [0], "k+", ms=7)
    ax.text(0.5 * k, -0.24, r"$\mathbf{V}$", ha="center", fontsize=11)

    theta0 = np.radians(115.0)  # 静止系での向きの例
    tip = np.array([k + np.cos(theta0), np.sin(theta0)])
    ax.annotate(
        "",
        xy=tip,
        xytext=(k, 0),
        arrowprops={"arrowstyle": "->", "color": "tab:orange", "lw": 1.6},
    )
    ax.annotate(
        "", xy=tip, xytext=(0, 0), arrowprops={"arrowstyle": "->", "color": "tab:blue", "lw": 1.8}
    )
    ax.text(*(0.5 * (tip + [k, 0]) + [0.08, 0.0]), r"$\mathbf{v}_0$", color="tab:orange")
    ax.text(*(0.55 * tip + [-0.22, 0.04]), r"$\mathbf{v}$", color="tab:blue")

    if k > 1.0:
        limit = max_lab_angle(k)
        reach = k + 1.15
        for sign in (1.0, -1.0):
            ax.plot(
                [0, reach * np.cos(limit)],
                [0, sign * reach * np.sin(limit)],
                ls="--",
                color="tab:red",
                lw=1,
            )
        ax.text(0.45, -0.55, r"$\theta_{\max}$", color="tab:red", fontsize=11)
        # 同じ向き θ に 2 つの速度がある
        theta = 0.6 * limit
        root = np.sqrt(1.0 - (k * np.sin(theta)) ** 2)
        for speed in (k * np.cos(theta) + root, k * np.cos(theta) - root):
            ax.plot(speed * np.cos(theta), -speed * np.sin(theta), "o", color="tab:green", ms=6)
        far = k * np.cos(theta) + root
        ax.plot([0, far * np.cos(theta)], [0, -far * np.sin(theta)], color="tab:green", lw=0.9)
    ax.set_aspect("equal")
    ax.set_xlim(-1.25, 3.4)
    ax.set_ylim(-1.35, 1.35)
    ax.set_xlabel(r"$v_x/v_0$")


def render_velocity_figure(ratios: Sequence[float], path: Path) -> Path:
    """親が遅い場合と速い場合の、破片の速度の取り得る範囲。"""
    setup_matplotlib()
    fig, axes = plt.subplots(1, len(ratios), figsize=(6.4, 2.6), constrained_layout=True)
    for ax, k in zip(axes, ratios, strict=True):
        _velocity_panel(ax, k)
        title = "どの向きにも飛ぶ" if k < 1 else "円すいの中にだけ飛ぶ"
        ax.set_title(f"$V/v_0={k:g}$：{title}", fontsize=10)
    axes[0].set_ylabel(r"$v_y/v_0$")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_distribution_figure(
    cases: Sequence[tuple[float, Histogram, Histogram]], path: Path
) -> Path:
    """左: 実験室系の角の分布。右: 運動エネルギーの分布。点は乱数実験、線は理論。

    ``cases`` は ``(V/v0, 角の分布, エネルギーの分布)``。
    """
    setup_matplotlib()
    fig, (ax_angle, ax_energy) = plt.subplots(1, 2, figsize=(6.4, 3.1), constrained_layout=True)
    for (k, angle, energy), color in zip(cases, COLORS, strict=False):
        label = f"$V/v_0={k:g}$"
        scale = np.radians(1.0) / angle.widths  # 1 度あたりの確率
        ax_angle.errorbar(
            np.degrees(angle.centers),
            angle.measured * scale,
            yerr=angle.error * scale,
            fmt=".",
            color=color,
            ms=4,
            lw=0.8,
        )
        fine = np.linspace(0.0, max_lab_angle(k) * (1 - 1e-4), 600)
        ax_angle.plot(
            np.degrees(fine),
            angle_density(fine, k) * np.radians(1.0),
            color=color,
            lw=1,
            label=label,
        )

        density = energy.measured / energy.widths
        ax_energy.errorbar(
            energy.centers,
            density,
            yerr=energy.error / energy.widths,
            fmt=".",
            color=color,
            ms=4,
            lw=0.8,
        )
        ax_energy.hlines(
            1.0 / (energy.edges[-1] - energy.edges[0]),
            energy.edges[0],
            energy.edges[-1],
            color=color,
            lw=1,
            label=label,
        )
    ax_angle.set_xlim(0, 180)
    ax_angle.set_xticks(range(0, 181, 45))
    ax_angle.set_yscale("log")
    ax_angle.set_ylim(2e-4, 1.0)
    ax_angle.set_xlabel(r"実験室系の角 $\theta$ [度]")
    ax_angle.set_ylabel("確率密度 [1/度]")
    ax_angle.set_title("向き：親が速いと前に集まる", fontsize=10)
    ax_angle.legend(fontsize=8, loc="upper right")
    ax_energy.set_xlabel("破片の運動エネルギー")
    ax_energy.set_ylabel("確率密度")
    ax_energy.set_title("エネルギー：一様に分布する", fontsize=10)
    ax_energy.set_ylim(0, 0.045)
    ax_energy.legend(fontsize=8, loc="upper right")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_burst_video(bursts: Sequence[Burst], path: Path, fps: int = 30) -> Path:
    """親が動きながら崩壊し、破片が円（球面の断面）をつくって広がる（実時間）。"""
    setup_matplotlib()
    fig, axes = plt.subplots(len(bursts), 1, figsize=(8.0, 3.6 * len(bursts)))
    axes = np.atleast_1d(axes)
    t = bursts[0].t
    phi = np.linspace(0.0, 2.0 * np.pi, 200)
    artists = []
    for ax, burst in zip(axes, bursts, strict=True):
        k = burst.v_parent / burst.v0
        ax.set_aspect("equal")
        ax.set_xlim(-3.4, 10.6)
        ax.set_ylim(-3.4, 3.4)
        ax.set_title(f"親の速さ / 静止系での破片の速さ = {k:g}")
        ax.plot([0], [0], "kx", ms=8)
        if k > 1.0:
            limit = max_lab_angle(k)
            for sign in (1.0, -1.0):
                ax.plot(
                    [0, 12 * np.cos(limit)],
                    [0, sign * 12 * np.sin(limit)],
                    ls="--",
                    color="tab:red",
                    lw=1,
                )
        (envelope,) = ax.plot([], [], color="0.7", lw=1)
        (dots,) = ax.plot([], [], ".", color="tab:blue", ms=6)
        (parent,) = ax.plot([], [], "o", color="k", ms=10)
        artists.append((burst, envelope, dots, parent))
    clock = axes[0].text(0.02, 0.9, "", transform=axes[0].transAxes)
    fig.tight_layout()
    dt = float(t[1] - t[0])
    frames = np.round((np.arange(t[0], t[-1] + 1e-9, 1.0 / fps) - t[0]) / dt).astype(int)

    def update(i: int):
        for burst, envelope, dots, parent in artists:
            now = t[i]
            if now >= 0:
                envelope.set_data(
                    burst.parent_x[i] + burst.v0 * now * np.cos(phi),
                    burst.v0 * now * np.sin(phi),
                )
                dots.set_data(burst.fragments[i, :, 0], burst.fragments[i, :, 1])
                parent.set_data([], [])
            else:
                envelope.set_data([], [])
                dots.set_data([], [])
                parent.set_data([burst.parent_x[i]], [0.0])
        clock.set_text(f"t = {t[i]:+.2f} s（0 で崩壊）")
        return ()

    anim = animation.FuncAnimation(fig, update, frames=frames, blit=False)
    saved = save_animation(anim, path, fps=fps)
    plt.close(fig)
    return saved
