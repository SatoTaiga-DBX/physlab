"""ガリレイの相対性の図と動画（電車の中の投げ上げ、作用のずれ、二体問題の 2 つの慣性系）。

設計書: docs/design/005-analytical-mechanics-ch1.md（§3.3, §3.4, §5）
結果データ（application）を受け取って描くだけで、計算はしない。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

from physlab.application.mechanics.galilean import ActionShift, TrainThrow, TwoBodyFrames
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_animation, save_figure

HEAVY = "tab:blue"  # 重い方の天体
LIGHT = "tab:orange"
ACCENT = "#c0392b"


def render_train_figure(r: TrainThrow, path: Path, n_marks: int = 9) -> Path:
    """電車の中で球を真上に投げる。電車の人には直線、地面の人には放物線。加速度はどちらも同じ。"""
    setup_matplotlib()
    fig, (ax_k, ax_g) = plt.subplots(
        1, 2, figsize=(6.4, 2.9), gridspec_kw={"width_ratios": [1, 2.3]}, constrained_layout=True
    )
    marks = np.linspace(0, r.t.size - 1, n_marks).astype(int)
    apex = float(r.train[:, 1].max())
    car_w, car_h = 3.0, 0.9  # 電車の窓枠の大きさ

    for ax, pos, shift_car in ((ax_k, r.train, 0.0), (ax_g, r.ground, r.train_speed)):
        ax.plot(pos[:, 0], pos[:, 1], "-", color="0.6", lw=1)
        ax.scatter(
            pos[marks, 0],
            pos[marks, 1],
            c=r.t[marks],
            cmap="viridis",
            s=34,
            zorder=5,
            edgecolor="k",
            linewidth=0.5,
        )
        # 電車（始めと終わり）。電車系では動かない
        for t_car in (0.0, r.t[-1]):
            x0 = shift_car * t_car - car_w / 2
            ax.add_patch(Rectangle((x0, -car_h), car_w, car_h, fc="0.88", ec="0.4", lw=1, zorder=1))
        ax.set_ylim(-1.2 * car_h, apex + 0.9)
        ax.set_ylabel("高さ $y$ [m]")
    ax_k.set_xlim(-2.2, 2.2)
    ax_g.set_xlim(-2.2, r.ground[-1, 0] + 2.2)
    ax_k.set_xlabel("$x'$ [m]")
    ax_g.set_xlabel("$x$ [m]")
    ax_k.set_title("電車の中の人（$K'$ 系）", fontsize=10)
    ax_g.set_title("地面の人（$K$ 系）", fontsize=10)
    ax_k.text(
        0.5,
        0.97,
        "まっすぐ上へ、そして落ちる",
        transform=ax_k.transAxes,
        ha="center",
        va="top",
        fontsize=8,
    )
    ax_g.text(
        0.5,
        0.97,
        "放物線。電車は右へ $V$ で動く",
        transform=ax_g.transAxes,
        ha="center",
        va="top",
        fontsize=8,
    )
    ax_g.annotate(
        "",
        (r.ground[-1, 0] + 1.2, -car_h / 2),
        (r.ground[-1, 0] - 0.3, -car_h / 2),
        arrowprops={"arrowstyle": "->", "color": "0.3"},
    )
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_action_shift_figure(s: ActionShift, path: Path) -> Path:
    """両端を固定した経路の族で、2 つの慣性系の作用は定数だけずれる。"""
    setup_matplotlib()
    fig, (ax_s, ax_d) = plt.subplots(1, 2, figsize=(6.4, 2.8), constrained_layout=True)
    ax_s.plot(s.eps, s.action_k, color="tab:blue", lw=2, label="$S$（$K$ 系）")
    ax_s.plot(s.eps, s.action_kprime, color="tab:green", lw=2, label="$S'$（$K'$ 系）")
    mid = int(np.argmin(s.action_k))
    ax_s.plot([s.eps[mid]], [s.action_k[mid]], "*", color="tab:blue", ms=12)
    ax_s.plot([s.eps[mid]], [s.action_kprime[mid]], "*", color="tab:green", ms=12)
    ax_s.set_xlabel("ずらし幅 $\\varepsilon$ [m]")
    ax_s.set_ylabel("作用 [J·s]")
    ax_s.set_title("曲線は平行にずれるだけ（最小の位置は同じ）", fontsize=9.5)
    ax_s.legend(fontsize=8.5)

    ax_d.plot(s.eps, s.shift, "o-", color="tab:purple", ms=4, label="$S-S'$（数値）")
    ax_d.axhline(
        s.shift_theory,
        color=ACCENT,
        ls="--",
        lw=1.2,
        label=f"$mV\\Delta q'+\\frac{{1}}{{2}}mV^2\\tau={s.shift_theory:.1f}$",
    )
    centre = s.shift_theory
    ax_d.set_ylim(centre - 6, centre + 6)
    ax_d.set_xlabel("ずらし幅 $\\varepsilon$ [m]")
    ax_d.set_ylabel("$S-S'$ [J·s]")
    ax_d.set_title("差は経路によらない定数", fontsize=9.5)
    ax_d.legend(fontsize=8.5, loc="upper right")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def _body_sizes(frames: TwoBodyFrames) -> tuple[float, float]:
    """質量に応じた点の大きさ（面積が質量に比例）。"""
    return 16.0, 16.0 * (frames.params.m2 / frames.params.m1) ** (1 / 3)


def _draw_two_body_panel(ax, q: np.ndarray, centre: np.ndarray | None, title: str) -> None:
    ax.plot(q[:, 0], q[:, 1], color=HEAVY, lw=1.1, alpha=0.9)
    ax.plot(q[:, 2], q[:, 3], color=LIGHT, lw=1.1, alpha=0.9)
    if centre is not None:
        ax.plot(centre[:, 0], centre[:, 1], ":", color="k", lw=1.2)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")


def render_two_body_figure(frames: TwoBodyFrames, path: Path) -> Path:
    """重心系 ``K'`` は重心が静止、実験室系 ``K`` は重心が等速で動く。同じ運動の 2 つの見え方。"""
    setup_matplotlib()
    fig, (ax_cm, ax_lab) = plt.subplots(
        1, 2, figsize=(6.4, 2.8), gridspec_kw={"width_ratios": [1, 1.9]}, constrained_layout=True
    )
    _draw_two_body_panel(ax_cm, frames.q_cm, None, "重心系 $K'$：重心は止まっている")
    ax_cm.plot([0], [0], "k+", ms=9)
    _draw_two_body_panel(
        ax_lab, frames.q_lab, frames.center_lab, "実験室系 $K$：重心が動く（点線）"
    )

    # 同じ時刻の姿勢を、両方の系に重ねて描く（二つの錘と、それを結ぶ線）
    for frac in (0.0, 0.25, 0.5, 0.75):
        i = int(frac * (frames.t.size - 1))
        for ax, q in ((ax_cm, frames.q_cm), (ax_lab, frames.q_lab)):
            ax.plot([q[i, 0], q[i, 2]], [q[i, 1], q[i, 3]], "-", color="0.45", lw=0.8)
            ax.plot([q[i, 0]], [q[i, 1]], "o", color=HEAVY, ms=6, mec="k", mew=0.5, zorder=5)
            ax.plot([q[i, 2]], [q[i, 3]], "o", color=LIGHT, ms=4.5, mec="k", mew=0.5, zorder=5)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_two_body_video(frames: TwoBodyFrames, path: Path, fps: int = 30) -> Path:
    """二体運動を重心系 ``K'``（左）と実験室系 ``K``（右）で並べた動画。再生は実時間と一致する。"""
    setup_matplotlib()
    dt = float(frames.t[1] - frames.t[0])
    idx = np.round(np.arange(0.0, frames.t[-1], 1.0 / fps) / dt).astype(int)
    size1, size2 = _body_sizes(frames)

    fig, (ax_cm, ax_lab) = plt.subplots(
        1, 2, figsize=(10, 4.6), gridspec_kw={"width_ratios": [1, 1.9]}
    )
    for ax, q, title in (
        (ax_cm, frames.q_cm, "重心系 $K'$（重心は止まっている）"),
        (ax_lab, frames.q_lab, "実験室系 $K$（重心が等速で動く）"),
    ):
        xs, ys = np.concatenate([q[:, 0], q[:, 2]]), np.concatenate([q[:, 1], q[:, 3]])
        ax.set_xlim(xs.min() - 0.4, xs.max() + 0.4)
        ax.set_ylim(ys.min() - 0.4, ys.max() + 0.4)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.grid(False)
    ax_cm.plot([0], [0], "k+", ms=10)
    ax_lab.plot(
        frames.center_lab[:, 0], frames.center_lab[:, 1], ":", color="k", lw=1, label="重心"
    )
    ax_lab.legend(loc="upper left")

    def artists(ax):
        (tr1,) = ax.plot([], [], "-", color=HEAVY, lw=1, alpha=0.7)
        (tr2,) = ax.plot([], [], "-", color=LIGHT, lw=1, alpha=0.7)
        (link,) = ax.plot([], [], "-", color="0.5", lw=0.8)
        (b1,) = ax.plot([], [], "o", color=HEAVY, ms=size1, mec="k")
        (b2,) = ax.plot([], [], "o", color=LIGHT, ms=size2, mec="k")
        return tr1, tr2, link, b1, b2

    cm_art, lab_art = artists(ax_cm), artists(ax_lab)
    time_text = ax_cm.text(0.03, 0.04, "", transform=ax_cm.transAxes)
    fig.tight_layout()

    def update(i: int):
        out = []
        for art, q in ((cm_art, frames.q_cm), (lab_art, frames.q_lab)):
            tr1, tr2, link, b1, b2 = art
            tr1.set_data(q[: i + 1, 0], q[: i + 1, 1])
            tr2.set_data(q[: i + 1, 2], q[: i + 1, 3])
            link.set_data([q[i, 0], q[i, 2]], [q[i, 1], q[i, 3]])
            b1.set_data([q[i, 0]], [q[i, 1]])
            b2.set_data([q[i, 2]], [q[i, 3]])
            out += [tr1, tr2, link, b1, b2]
        time_text.set_text(f"t = {frames.t[i]:5.2f}")
        return (*out, time_text)

    anim = animation.FuncAnimation(fig, update, frames=idx, blit=True)
    saved = save_animation(anim, path, fps=fps)
    plt.close(fig)
    return saved
