"""最小作用の原理の図と動画。結果データ（application）を受け取って描くだけで、計算はしない。

設計書: docs/design/005-analytical-mechanics-ch1.md（§3.2, §5）
教科書に貼る図は PDF（ベクタ）で、幅は本文幅（約 6.3 インチ）に合わせる。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import animation, colors
from matplotlib import pyplot as plt

from physlab.application.mechanics.least_action import (
    ActionLandscape,
    ActionScan,
    CurvatureScan,
    PathRelaxation,
)
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_animation, save_figure

TRUE_COLOR = "#16325c"  # 真の運動（停留経路）
ACCENT = "#c0392b"


def render_path_family_figure(scan: ActionScan, path: Path) -> Path:
    """左: 両端を固定した経路の族。右: それぞれの作用 ``S(ε)``（真の運動で最小）。"""
    setup_matplotlib()
    fig, (ax_p, ax_s) = plt.subplots(1, 2, figsize=(6.4, 2.9), constrained_layout=True)
    t1, t2 = scan.t[0], scan.t[-1]
    cmap = plt.get_cmap("coolwarm")
    shown = np.linspace(0, scan.eps.size - 1, 5).astype(int)
    norm = colors.Normalize(scan.eps.min(), scan.eps.max())

    for i in shown:
        if scan.eps[i] != 0:
            ax_p.plot(scan.t, scan.paths[i], color=cmap(norm(scan.eps[i])), lw=1.2)
    ax_p.plot(scan.t, scan.true_path, color=TRUE_COLOR, lw=2.6, label="真の運動 ($\\varepsilon=0$)")
    ax_p.plot([t1, t2], [scan.true_path[0], scan.true_path[-1]], "ko", ms=5, zorder=5)
    ax_p.set_xlabel("時間 $t$ [s]")
    ax_p.set_ylabel("高さ $q$ [m]")
    ax_p.set_title("両端を固定した経路の族", fontsize=10)
    ax_p.legend(loc="lower center", fontsize=8)

    fine = np.linspace(scan.eps.min(), scan.eps.max(), 200)
    ax_s.plot(
        fine, scan.action_true + scan.fit_curvature * fine**2 + scan.fit_linear * fine, color="0.35"
    )
    ax_s.scatter(
        scan.eps[shown],
        scan.action[shown],
        c=[cmap(norm(e)) for e in scan.eps[shown]],
        s=40,
        zorder=5,
        edgecolor="k",
        linewidth=0.6,
    )
    ax_s.scatter([0], [scan.action_true], s=70, marker="*", color=TRUE_COLOR, zorder=6)
    ax_s.annotate(
        f"$S_*={scan.action_true:.2f}$",
        (0, scan.action_true),
        (0.15, 0.12),
        textcoords="axes fraction",
        fontsize=9,
        arrowprops={"arrowstyle": "-", "color": "0.4"},
    )
    ax_s.set_xlabel("ずらし幅 $\\varepsilon$ [m]")
    ax_s.set_ylabel("作用 $S$ [J·s]")
    ax_s.set_title("作用 $S(\\varepsilon)$：真の運動で最小", fontsize=10)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_second_variation_figure(
    short: ActionScan, long_mode1: ActionScan, long_mode2: ActionScan, path: Path
) -> Path:
    """同じ調和振動子でも、時間が短ければ最小、長ければ鞍点（モードで上がる/下がる）。"""
    setup_matplotlib()
    fig, axes = plt.subplots(1, 3, figsize=(6.4, 2.5), constrained_layout=True)
    panels = [
        (short, f"短い区間 $\\tau={short.t[-1]:.1f}$ s\n第 1 モード：上がる"),
        (long_mode1, f"長い区間 $\\tau={long_mode1.t[-1]:.1f}$ s\n第 1 モード：下がる"),
        (long_mode2, f"長い区間 $\\tau={long_mode2.t[-1]:.1f}$ s\n第 2 モード：上がる"),
    ]
    for ax, (scan, title) in zip(axes, panels, strict=True):
        delta = scan.action - scan.action_true
        ax.plot(scan.eps, delta, color=ACCENT if scan.theory_curvature < 0 else "tab:blue", lw=2)
        ax.axhline(0, color="0.5", lw=0.8)
        ax.plot([0], [0], "*", color=TRUE_COLOR, ms=11, zorder=5)
        ax.set_title(title, fontsize=8.5)
        ax.set_xlabel("$\\varepsilon$ [m]")
        # 係数の文字は、曲線のない側（上向きの碗なら上、下向きなら下）に置く
        ax.text(
            0.5,
            0.86 if scan.theory_curvature > 0 else 0.07,
            f"$c_{scan.mode}={scan.theory_curvature:+.2f}$",
            transform=ax.transAxes,
            ha="center",
            fontsize=9,
        )
    axes[0].set_ylabel("$S(\\varepsilon)-S_*$ [J·s]")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_curvature_figure(scan: CurvatureScan, path: Path) -> Path:
    """作用の曲率 ``c_n(τ)``。零点は ``τ = nπ/ω``（そこから先でモード n が最小でなくなる）。"""
    setup_matplotlib()
    fig, ax = plt.subplots(figsize=(6.0, 3.0), constrained_layout=True)
    omega = scan.system.omega
    palette = ["tab:blue", "tab:orange", "tab:green"]
    for color, n in zip(palette, scan.modes, strict=False):
        ax.plot(
            scan.tau, scan.theory[n], color=color, lw=1.8, label=f"第 {n} モード $c_{n}(\\tau)$"
        )
        ax.plot(scan.tau[::8], scan.numeric[n][::8], ".", color=color, ms=5)  # 数値積分の結果
        ax.axvline(n * np.pi / omega, color=color, ls=":", lw=1)
        zero = f"{n}\\pi/\\omega" if n > 1 else "\\pi/\\omega"
        ax.text(n * np.pi / omega + 0.01, -23.5, f"$\\tau={zero}$", fontsize=8.5, color=color)
    ax.axhline(0, color="0.3", lw=0.8)
    ax.axvspan(np.pi / omega, scan.tau.max(), color=ACCENT, alpha=0.06)
    ax.text(
        np.pi / omega + 0.02,
        53,
        "$\\tau>\\pi/\\omega$：第 1 モードの曲率が負（鞍点）",
        fontsize=8.5,
        color=ACCENT,
    )
    ax.set_ylim(-25, 60)
    ax.set_xlabel("時間の長さ $\\tau$ [s]")
    ax.set_ylabel("作用の曲率 $c_n$ [J·s/m²]")
    ax.set_title("曲率が負になると、真の運動は「最小」ではなくなる（点: 数値積分）", fontsize=10)
    ax.legend(loc="upper right", fontsize=8.5, bbox_to_anchor=(1.0, 0.9))
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_landscape_figure(bowl: ActionLandscape, saddle: ActionLandscape, path: Path) -> Path:
    """経路の空間（2 つのモードの振幅 a, b）の上の地形 ``S − S*``。碗と鞍。"""
    setup_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.0), constrained_layout=True)
    lowest = float(saddle.delta_action.min())
    saddle_norm = colors.TwoSlopeNorm(vcenter=0.0, vmin=lowest, vmax=-lowest)  # 0 を白にそろえる
    for ax, land, title, cmap, norm in (
        (axes[0], bowl, f"短い区間 $\\tau={bowl.tau:.1f}$ s：碗（最小）", "viridis", None),
        (axes[1], saddle, f"長い区間 $\\tau={saddle.tau:.1f}$ s：鞍", "RdBu_r", saddle_norm),
    ):
        mesh = ax.contourf(land.a, land.b, land.delta_action, levels=16, cmap=cmap, norm=norm)
        ax.contour(land.a, land.b, land.delta_action, levels=16, colors="k", linewidths=0.3)
        ax.plot([0], [0], "*", color="w", mec="k", ms=13, zorder=5)
        ax.set_aspect("equal")
        ax.set_xlabel("第 1 モードの振幅 $a$ [m]")
        ax.set_title(title, fontsize=9.5)
        ax.grid(False)
        fig.colorbar(mesh, ax=ax, shrink=0.8, label="$S-S_*$ [J·s]")
    axes[0].set_ylabel("第 2 モードの振幅 $b$ [m]")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def _action_excess(r: PathRelaxation) -> np.ndarray:
    """停留経路の作用との差 ``S_d − S_*``。対数軸に描くため、負にならないよう下限を置く。"""
    return np.maximum(r.action - r.action_true, 1e-12)


def render_relaxation_figure(r: PathRelaxation, path: Path) -> Path:
    """作用を下げていくと、でたらめな経路が真の運動（放物線）に収束する様子の静止画版。"""
    setup_matplotlib()
    fig, (ax_p, ax_s) = plt.subplots(1, 2, figsize=(6.4, 2.9), constrained_layout=True)
    wanted = [0, 5, 30, 200, 1000, int(r.iterations[-1])]
    picks = [int(np.argmin(np.abs(r.iterations - k))) for k in wanted]
    cmap = plt.get_cmap("viridis_r")
    for j, i in enumerate(dict.fromkeys(picks)):
        ax_p.plot(
            r.t,
            r.paths[i],
            color=cmap(j / (len(picks) - 1)),
            lw=1.3,
            label=f"反復 {int(r.iterations[i])}",
        )
    ax_p.plot(r.t, r.true_path, "--", color="k", lw=1.4, label="放物線（解析解）")
    ax_p.set_xlabel("時間 $t$ [s]")
    ax_p.set_ylabel("高さ $q$ [m]")
    ax_p.set_title("作用を下げ続けた経路", fontsize=10)
    low, high = ax_p.get_ylim()
    ax_p.set_ylim(low, high + 0.55 * (high - low))  # 凡例を置く余白
    ax_p.legend(fontsize=6.5, loc="upper center", ncol=3)

    excess = _action_excess(r)
    ax_s.plot(r.iterations, excess, color=TRUE_COLOR, lw=1.8)
    ax_s.set_xscale("symlog", linthresh=1)
    ax_s.set_yscale("log")
    ax_s.set_xlim(0, r.iterations[-1] * 1.5)
    ax_s.set_xlabel("反復回数")
    ax_s.set_ylabel("作用の超過分 $S_d-S_*$ [J·s]")
    ax_s.set_title("作用は単調に下がる", fontsize=10)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_relaxation_video(
    r: PathRelaxation, path: Path, fps: int = 24, repeat: int = 2, hold_seconds: float = 1.5
) -> Path:
    """経路が反復とともに真の運動へ収束していく動画。各フレームを ``repeat`` 回ずつ見せる。"""
    setup_matplotlib()
    fig, (ax_p, ax_s) = plt.subplots(
        1, 2, figsize=(10, 4.8), gridspec_kw={"width_ratios": [1.25, 1]}
    )
    pad = 0.08 * np.ptp(r.paths)
    ax_p.set_xlim(r.t[0], r.t[-1])
    ax_p.set_ylim(r.paths.min() - pad, r.paths.max() + pad)
    ax_p.set_xlabel("時間 $t$ [s]")
    ax_p.set_ylabel("高さ $q$ [m]")
    ax_p.set_title("経路 $q(t)$（両端は固定）")
    ax_p.plot(r.t, r.true_path, "--", color="0.45", lw=1.6, label="真の運動（放物線）")
    (line,) = ax_p.plot([], [], "-o", color=TRUE_COLOR, lw=2, ms=3.5, label="いまの経路")
    ax_p.plot([r.t[0], r.t[-1]], [r.paths[0, 0], r.paths[0, -1]], "ko", ms=7, zorder=5)
    ax_p.legend(loc="lower right")  # 経路が通らない右下に置く
    text = ax_p.text(0.03, 0.04, "", transform=ax_p.transAxes, fontsize=11)

    excess = _action_excess(r)
    ax_s.plot(r.iterations, excess, color="0.8", lw=1)
    ax_s.set_xscale("symlog", linthresh=1)
    ax_s.set_yscale("log")
    ax_s.set_xlim(0, r.iterations[-1] * 1.5)
    ax_s.set_ylim(1e-10, 1e3)
    ax_s.set_xlabel("反復回数")
    ax_s.set_ylabel("作用の超過分 $S_d-S_*$ [J·s]")
    ax_s.set_title("作用が下がると、経路は真の運動に近づく")
    (trace,) = ax_s.plot([], [], color=TRUE_COLOR, lw=2)
    (dot,) = ax_s.plot([], [], "o", color=TRUE_COLOR, ms=8)
    fig.tight_layout()

    n = r.iterations.size
    order = np.repeat(np.arange(n), repeat).tolist() + [n - 1] * int(hold_seconds * fps)

    def update(i: int):
        line.set_data(r.t, r.paths[i])
        trace.set_data(r.iterations[: i + 1], excess[: i + 1])
        dot.set_data([r.iterations[i]], [excess[i]])
        text.set_text(
            f"反復 {int(r.iterations[i])}\n$S_d$ = {r.action[i]:.2f}   "
            f"真の運動との最大のずれ = {r.max_error[i]:.2g} m"
        )
        return line, trace, dot, text

    anim = animation.FuncAnimation(fig, update, frames=order, blit=True)
    saved = save_animation(anim, path, fps=fps)
    plt.close(fig)
    return saved
