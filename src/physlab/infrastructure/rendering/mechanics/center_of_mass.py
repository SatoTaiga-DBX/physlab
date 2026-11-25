"""慣性中心（重心）の図と動画。結果データを受け取って描くだけで、計算はしない。

設計書: docs/design/006-conservation-laws.md
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt

from physlab.application.mechanics.center_of_mass import ThrowResult
from physlab.domain.mechanics.cluster import as_vectors
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_animation, save_figure

PARTICLE_COLORS = ("tab:blue", "tab:orange", "tab:green")
LOOP = [0, 1, 2, 0]


def render_throw_figure(result: ThrowResult, path: Path, n_snapshots: int = 8) -> Path:
    """各質点は複雑に回るが、重心は放物線を描く。"""
    setup_matplotlib()
    r = as_vectors(result.q)
    fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    for a, color in enumerate(PARTICLE_COLORS):
        ax.plot(r[:, a, 0], r[:, a, 1], color=color, lw=1, alpha=0.8, label=f"質点 {a + 1}")
    for k, i in enumerate(np.linspace(0, result.t.size - 1, n_snapshots).astype(int)):
        ax.plot(r[i, LOOP, 0], r[i, LOOP, 1], "-", color="0.55", lw=0.8)
        for a, color in enumerate(PARTICLE_COLORS):
            ax.plot(r[i, a, 0], r[i, a, 1], "o", color=color, ms=4 + 1.5 * k / n_snapshots)
        ax.plot(*result.center[i], "k+", ms=9, mew=1.8, zorder=6)
    ax.plot(result.center[:, 0], result.center[:, 1], color="k", lw=2.2, label="重心")
    ax.plot(
        result.center_exact[:, 0],
        result.center_exact[:, 1],
        "--",
        color="gold",  # 黒い重心の線の上でも、凡例でも見える色
        lw=1.1,
        label="放物線の式",
        zorder=4,
    )
    ax.set_aspect("equal")
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")
    ax.set_title("回りながら振動して飛ぶ 3 質点（+ は同じ時刻の重心）", fontsize=10)
    ax.legend(fontsize=7.5, loc="upper left", ncol=2)
    low, high = ax.get_ylim()
    ax.set_ylim(low, high + 0.25 * (high - low))
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_decomposition_figure(result: ThrowResult, path: Path) -> Path:
    """エネルギーと角運動量を「重心の分」と「重心に対する分」に分ける。"""
    setup_matplotlib()
    fig, (ax_e, ax_l) = plt.subplots(1, 2, figsize=(6.4, 2.9), constrained_layout=True)
    t = result.t
    ax_e.plot(t, result.kinetic_center + result.gravity, color="k", lw=2, label="重心の運動 + 重力")
    ax_e.plot(t, result.internal_energy, color="tab:red", lw=2, label="内部エネルギー")
    ax_e.plot(t, result.kinetic_relative, color="tab:orange", lw=0.9, label="内部の運動")
    ax_e.plot(t, result.spring, color="tab:purple", lw=0.9, label="ばね")
    ax_e.set_ylim(
        0,
        1.55
        * max(
            result.internal_energy.max(), result.total_energy.max() - result.internal_energy.min()
        ),
    )
    ax_e.set_xlabel("時間 $t$ [s]")
    ax_e.set_ylabel("エネルギー [J]")
    ax_e.set_title("エネルギーの分け方", fontsize=10)
    ax_e.legend(fontsize=7, loc="upper center", ncol=2)

    ax_l.plot(t, result.angular_total, color="k", lw=2, label="原点まわり（全体）")
    ax_l.plot(t, result.angular_center, "--", color="0.45", lw=1.6, label=r"重心の分 $R\times P$")
    ax_l.plot(t, result.angular_relative, color="tab:red", lw=2, label="重心まわり")
    ax_l.set_xlabel("時間 $t$ [s]")
    ax_l.set_ylabel("角運動量 [kg·m²/s]")
    ax_l.set_title("重力があっても重心まわりは一定", fontsize=10)
    ax_l.legend(fontsize=7.5, loc="lower left")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_throw_video(result: ThrowResult, path: Path, fps: int = 30, slow: float = 4.0) -> Path:
    """投げた 3 質点の動画。飛行が 2 秒足らずなので、``1/slow`` 倍速で見せる（画面に表示する）。"""
    setup_matplotlib()
    r = as_vectors(result.q)
    dt = float(result.t[1] - result.t[0])
    idx = np.round(np.arange(0.0, result.t[-1], 1.0 / (fps * slow)) / dt).astype(int)

    fig, ax = plt.subplots(figsize=(9, 5.6))
    pad = 0.4
    ax.set_xlim(r[..., 0].min() - pad, r[..., 0].max() + pad)
    ax.set_ylim(r[..., 1].min() - pad, r[..., 1].max() + 1.2)
    ax.set_aspect("equal")
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")
    ax.set_title(f"重心は放物線を描く（1/{slow:g} 倍速）")
    ax.plot(
        result.center_exact[:, 0],
        result.center_exact[:, 1],
        ":",
        color="0.5",
        lw=1.2,
        label="放物線の式",
    )
    trails = [ax.plot([], [], color=c, lw=1, alpha=0.6)[0] for c in PARTICLE_COLORS]
    (triangle,) = ax.plot([], [], "-", color="0.5", lw=1.2)
    dots = [ax.plot([], [], "o", color=c, ms=11)[0] for c in PARTICLE_COLORS]
    (cm_trail,) = ax.plot([], [], color="k", lw=2, label="重心の軌跡")
    (cm_dot,) = ax.plot([], [], "k+", ms=13, mew=2.5)
    text = ax.text(0.02, 0.93, "", transform=ax.transAxes)
    ax.legend(loc="upper right")
    fig.tight_layout()

    def update(i: int):
        for a, (trail, dot) in enumerate(zip(trails, dots, strict=True)):
            trail.set_data(r[: i + 1, a, 0], r[: i + 1, a, 1])
            dot.set_data([r[i, a, 0]], [r[i, a, 1]])
        triangle.set_data(r[i, LOOP, 0], r[i, LOOP, 1])
        cm_trail.set_data(result.center[: i + 1, 0], result.center[: i + 1, 1])
        cm_dot.set_data([result.center[i, 0]], [result.center[i, 1]])
        text.set_text(f"t = {result.t[i]:.2f} s")
        return (*trails, triangle, *dots, cm_trail, cm_dot, text)

    anim = animation.FuncAnimation(fig, update, frames=idx, blit=True)
    saved = save_animation(anim, path, fps=fps)
    plt.close(fig)
    return saved
