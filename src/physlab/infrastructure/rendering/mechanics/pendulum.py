"""単振り子の動画・図の描画。結果データ（application）を受け取って描くだけで、計算はしない。

設計書: docs/design/002-pendulum.md（「可視化の設計」）
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt

from physlab.application.mechanics.pendulum import (
    EnergyComparison,
    PeriodCurve,
    PhasePortrait,
    SmallAngleComparison,
)
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_animation, save_figure

ACCENT = "#c0392b"  # 分離線の色
_METHOD_LABELS = {
    "euler": "前進オイラー（1次）",
    "rk4": "ルンゲ＝クッタ4次",
    "verlet": "速度ベルレ（2次・シンプレクティック）",
}


def render_comparison_video(result: SmallAngleComparison, path: Path, fps: int = 30) -> Path:
    """非線形解と小角近似解を並べた動画を書き出し、実際のパスを返す。"""
    t, theta, theta_lin = result.t, result.theta, result.theta_linear
    L = result.params.length
    theta0 = result.params.theta0
    duration = float(t[-1])
    dt = float(t[1] - t[0])
    # 動画の k コマ目の時刻 k/fps に最も近い計算ステップを使う（再生が実時間と一致する）
    frames = np.round(np.arange(0.0, duration, 1.0 / fps) / dt).astype(int)

    setup_matplotlib()
    fig, (ax_p, ax_t) = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={"width_ratios": [1, 1.4]})

    # 左: 振り子そのもの
    ax_p.set_xlim(-1.3 * L, 1.3 * L)
    # 90° を超えて振れると錘が支点より上に来るので、その分だけ上側の表示範囲を広げる
    y_top = max(0.3, 1.2 * -np.cos(theta0))
    ax_p.set_ylim(-1.3 * L, y_top * L)
    ax_p.set_aspect("equal")
    ax_p.set_title("振り子の動き")
    ax_p.grid(False)
    ax_p.plot([0], [0], "k^", ms=8)
    (rod_nl,) = ax_p.plot([], [], "-", color="tab:blue", lw=2, label="厳密（非線形）")
    (bob_nl,) = ax_p.plot([], [], "o", color="tab:blue", ms=16)
    (rod_lin,) = ax_p.plot([], [], "-", color="tab:orange", lw=2, alpha=0.6, label="小角近似")
    (bob_lin,) = ax_p.plot([], [], "o", color="tab:orange", ms=12, alpha=0.6)
    ax_p.legend(loc="upper right", fontsize=9)

    # 右: θ(t)
    ax_t.plot(t, np.rad2deg(theta), color="tab:blue", lw=1, alpha=0.35)
    ax_t.plot(t, np.rad2deg(theta_lin), color="tab:orange", lw=1, alpha=0.35)
    (dot_nl,) = ax_t.plot([], [], "o", color="tab:blue")
    (dot_lin,) = ax_t.plot([], [], "o", color="tab:orange")
    (trace_nl,) = ax_t.plot([], [], color="tab:blue", lw=2)
    (trace_lin,) = ax_t.plot([], [], color="tab:orange", lw=2, alpha=0.7)
    ax_t.set_xlim(0, duration)
    ax_t.set_ylim(-1.15 * np.rad2deg(abs(theta0)), 1.15 * np.rad2deg(abs(theta0)))
    ax_t.set_xlabel("時間 t [s]")
    ax_t.set_ylabel("角度 θ [deg]")
    ax_t.set_title(
        f"初期角 {np.rad2deg(theta0):.0f}°：周期は厳密 {result.period_exact:.3f} s"
        f" ／ 小角近似 {result.period_linear:.3f} s"
    )
    time_text = ax_p.text(0.03, 0.05, "", transform=ax_p.transAxes)
    fig.tight_layout()

    def update(i: int):
        for rod, bob, th in ((rod_nl, bob_nl, theta[i]), (rod_lin, bob_lin, theta_lin[i])):
            x, y = L * np.sin(th), -L * np.cos(th)
            rod.set_data([0, x], [0, y])
            bob.set_data([x], [y])
        dot_nl.set_data([t[i]], [np.rad2deg(theta[i])])
        dot_lin.set_data([t[i]], [np.rad2deg(theta_lin[i])])
        trace_nl.set_data(t[: i + 1], np.rad2deg(theta[: i + 1]))
        trace_lin.set_data(t[: i + 1], np.rad2deg(theta_lin[: i + 1]))
        time_text.set_text(f"t = {t[i]:5.2f} s")
        return rod_nl, bob_nl, rod_lin, bob_lin, dot_nl, dot_lin, trace_nl, trace_lin, time_text

    anim = animation.FuncAnimation(fig, update, frames=frames, blit=True)
    saved = save_animation(anim, path, fps=fps)
    plt.close(fig)
    return saved


def render_energy_drift_figure(result: EnergyComparison, path: Path) -> Path:
    """積分法ごとのエネルギー誤差の推移を 1 枚の図にして保存する。"""
    setup_matplotlib()
    fig, ax = plt.subplots(figsize=(8, 5))
    for method, rel_err in result.relative_error.items():
        ax.plot(
            result.t,
            np.maximum(np.abs(rel_err), 1e-16),
            lw=1,
            label=_METHOD_LABELS.get(method, method),
        )
    ax.set_yscale("log")
    ax.set_ylim(1e-10, 1e1)
    ax.set_xlabel("時間 t [s]")
    ax.set_ylabel("エネルギーの相対誤差 |E - E0| / E0")
    ax.set_title(f"積分法ごとのエネルギー保存（刻み dt = {result.dt} s）")
    ax.legend(loc="lower right")
    fig.tight_layout()
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_period_ratio_figure(curve: PeriodCurve, path: Path) -> Path:
    """周期の比 T/T0 を初期角の関数として描く。厳密値と、冪展開の 1 項・2 項の近似。"""
    setup_matplotlib()
    fig, ax = plt.subplots(figsize=(6.0, 3.2), constrained_layout=True)
    ax.plot(curve.theta0_deg, curve.ratio_exact, color="tab:blue", lw=2.2, label="厳密（楕円積分）")
    ax.plot(
        curve.theta0_deg,
        curve.ratio_first,
        "--",
        color="tab:orange",
        lw=1.5,
        label=r"$1+\theta_0^2/16$",
    )
    ax.plot(
        curve.theta0_deg,
        curve.ratio_second,
        "--",
        color="tab:green",
        lw=1.5,
        label=r"$1+\theta_0^2/16+11\theta_0^4/3072$",
    )
    offsets = {30: (-6, 14), 120: (8, -16)}  # 注記の位置（点からのずらし）
    for deg in (30, 120):
        i = int(np.argmin(np.abs(curve.theta0_deg - deg)))
        ax.plot([deg], [curve.ratio_exact[i]], "o", color="tab:blue", ms=6)
        ax.annotate(
            f"{deg}°: {curve.ratio_exact[i]:.4f}",
            (deg, curve.ratio_exact[i]),
            offsets[deg],
            textcoords="offset points",
            fontsize=8.5,
        )
    ax.axhline(1.0, color="0.4", lw=0.8)
    ax.set_ylim(0.9, 4.2)
    ax.set_xlim(0, 180)
    ax.set_xlabel(r"初期角 $\theta_0$ [deg]")
    ax.set_ylabel(r"周期の比 $T/T_0$")
    ax.set_title("振幅が大きいほど周期は伸び、真上（180 度）で発散する", fontsize=10)
    ax.legend(loc="upper left", fontsize=8.5)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_phase_portrait_figure(portrait: PhasePortrait, path: Path) -> Path:
    """位相平面。エネルギー保存より、運動は E/m = 一定の等高線の上を動く。赤が分離線。"""
    setup_matplotlib()
    p = portrait.params
    fig, ax = plt.subplots(figsize=(6.2, 3.4), constrained_layout=True)
    levels = p.g * p.length * np.array([0.1, 0.4, 0.8, 1.2, 1.6, 1.9])  # 振動
    rotating = p.g * p.length * np.array([2.3, 3.0, 4.0, 5.0])  # 回転
    th = portrait.theta
    ax.contour(
        th, portrait.omega, portrait.energy, levels=levels, colors="tab:blue", linewidths=1.3
    )
    ax.contour(
        th, portrait.omega, portrait.energy, levels=rotating, colors="tab:green", linewidths=1.0
    )
    ax.contour(
        th,
        portrait.omega,
        portrait.energy,
        levels=[portrait.separatrix_energy],
        colors=ACCENT,
        linewidths=2.0,
    )
    label_box = {
        "fc": "white",
        "ec": "none",
        "alpha": 0.9,
        "pad": 1.5,
    }  # 等高線の上でも読めるように
    ax.text(0, 0, "振動", ha="center", va="center", fontsize=9, color="tab:blue", bbox=label_box)
    ax.text(0.5 * np.pi, 5.5, "回転", fontsize=9, color="tab:green", ha="center", bbox=label_box)
    ax.text(
        np.pi, 1.2, "分離線\n（真上で静止）", fontsize=8, color=ACCENT, ha="center", bbox=label_box
    )
    ax.set_xticks(np.pi * np.arange(-2, 3))
    ax.set_xticklabels([r"$-2\pi$", r"$-\pi$", "$0$", r"$\pi$", r"$2\pi$"])
    ax.set_xlabel(r"角度 $\theta$ [rad]")
    ax.set_ylabel(r"角速度 $\dot\theta$ [rad/s]")
    ax.set_title("単振り子の位相平面（$E/m$ が一定の曲線）", fontsize=10)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved
