"""単振り子の動画・図の描画。結果データ（application）を受け取って描くだけで、計算はしない。

設計書: docs/design/002-pendulum.md（「可視化の設計」）
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt

from physlab.application.mechanics.pendulum import EnergyComparison, SmallAngleComparison
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_animation, save_figure

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
