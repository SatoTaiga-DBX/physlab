"""弾性衝突の図と動画。結果データを受け取って描くだけで、計算はしない。

設計書: docs/design/007-collisions.md
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
from matplotlib.patches import Circle

from physlab.application.mechanics.elastic_collision import CollisionBatch, CollisionOutcome
from physlab.domain.mechanics.collision import (
    energy_transfer_fraction,
    lab_angles,
    momentum_circle,
)
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_animation, save_figure

NAMES = {
    "soft-disk": "柔らかい円板",
    "gaussian": "ガウス型",
    "yukawa": "湯川型",
    "steep-core": r"急な芯 $r^{-12}$",
}
STYLES = {
    "soft-disk": ("o", "tab:blue"),
    "gaussian": ("s", "tab:orange"),
    "yukawa": ("^", "tab:green"),
    "steep-core": ("D", "tab:red"),
}


def render_momentum_figure(
    cases: Sequence[tuple[float, float, Mapping[str, CollisionOutcome]]], path: Path
) -> Path:
    """質量比ごとに、衝突後の運動量を運動量の平面に打つ。``cases`` は ``(m1, m2, 結果)``。"""
    setup_matplotlib()
    fig, axes = plt.subplots(1, len(cases), figsize=(6.5, 2.9), constrained_layout=True)
    phi = np.linspace(0.0, 2.0 * np.pi, 300)
    for ax, (m1, m2, results) in zip(axes, cases, strict=True):
        speed = 1.0
        total = m1 * speed
        center, radius = momentum_circle(m1, m2, speed)
        # すべて入射の運動量 m1 v を単位にして描く
        # p2' が乗る円（実線）と、p1' = m1 v − p2' が乗る円（破線）
        for middle, style in ((center, "-"), (total - center, ":")):
            ax.plot(
                (middle + radius * np.cos(phi)) / total,
                radius * np.sin(phi) / total,
                style,
                color="0.55",
                lw=1,
            )
        ax.annotate(
            "", xy=(1.0, 0.0), xytext=(0.0, 0.0), arrowprops={"arrowstyle": "->", "color": "k"}
        )
        for name, result in results.items():
            marker, color = STYLES[name]
            p2 = result.p2 / total
            p1 = result.p1 / total
            ax.plot(p2[:, 0], p2[:, 1], marker, color=color, ms=4, ls="none", label=NAMES[name])
            ax.plot(p1[:, 0], p1[:, 1], marker, color=color, ms=4, mfc="none", ls="none")
        ax.set_aspect("equal")
        ax.set_xlim(-0.45, 1.55)
        ax.set_ylim(-1.0, 1.0)
        ax.set_title(f"$m_2/m_1={m2 / m1:g}$", fontsize=10)
        ax.text(-0.38, 0.84, "上: 入射（白抜き）", fontsize=7)
        ax.text(-0.38, -0.92, "下: 標的（塗り）", fontsize=7)
        ax.set_xlabel(r"$p_x/m_1v$")
    axes[0].set_ylabel(r"$p_y/m_1v$")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=4, fontsize=8)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def _center_of_mass(batch: CollisionBatch) -> np.ndarray:
    q = batch.q
    total = batch.m1 + batch.m2
    return (batch.m1 * q[..., 0:2] + batch.m2 * q[..., 2:4]) / total


def render_trajectory_figure(batches: Mapping[str, CollisionBatch], path: Path) -> Path:
    """同じ衝突パラメータで、力の形を変えた軌道。左: 実験室系、右: 重心系。"""
    setup_matplotlib()
    fig, (ax_lab, ax_cm) = plt.subplots(1, 2, figsize=(6.5, 3.2), constrained_layout=True)
    for name, batch in batches.items():
        _, color = STYLES[name]
        q = batch.q[:, 0]
        cm = _center_of_mass(batch)[:, 0]
        for ax, shift in ((ax_lab, np.zeros_like(cm)), (ax_cm, cm)):
            ax.plot(*(q[:, 0:2] - shift).T, color=color, lw=1.3, label=NAMES[name])
            ax.plot(*(q[:, 2:4] - shift).T, color=color, lw=1.3, ls="--")
    for ax in (ax_lab, ax_cm):
        ax.set_aspect("equal")
        ax.set_xlim(-4.5, 4.5)
        ax.set_ylim(-4.5, 4.5)
        ax.set_xlabel("$x$")
    ax_lab.set_ylabel("$y$")
    ax_lab.set_title("実験室系（実線: 入射、破線: 標的）", fontsize=9)
    ax_cm.set_title("重心系：向きが変わるだけ", fontsize=9)
    ax_lab.legend(fontsize=7, loc="lower left")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_mass_ratio_figure(points: Sequence[tuple[float, CollisionOutcome]], path: Path) -> Path:
    """左: 重心系の角と実験室系の角。右: 正面衝突で渡せるエネルギーの割合。

    ``points`` は ``(m2/m1, 数値解の結果)``（``m1 = 1``）。
    """
    setup_matplotlib()
    fig, (ax_angle, ax_energy) = plt.subplots(1, 2, figsize=(6.5, 3.0), constrained_layout=True)
    chi = np.linspace(0.0, np.pi, 400)
    for (ratio, result), color in zip(points, ("tab:red", "tab:blue", "tab:green"), strict=True):
        theta1, _ = lab_angles(chi, 1.0, ratio)
        ax_angle.plot(
            np.degrees(chi), np.degrees(theta1), color=color, lw=1.2, label=f"$m_2/m_1={ratio:g}$"
        )
        ax_angle.plot(np.degrees(result.chi), np.degrees(result.theta1), "o", color=color, ms=4)
    ax_angle.set_xlim(0, 180)
    ax_angle.set_ylim(0, 180)
    ax_angle.set_xticks(range(0, 181, 45))
    ax_angle.set_yticks(range(0, 181, 45))
    ax_angle.set_xlabel(r"重心系の角 $\chi$ [度]")
    ax_angle.set_ylabel(r"入射粒子の角 $\theta_1$ [度]")
    ax_angle.set_title("線: 公式、点: 運動方程式", fontsize=9)
    ax_angle.legend(fontsize=7.5, loc="upper left")

    ratio = np.geomspace(0.01, 300.0, 400)
    ax_energy.plot(ratio, energy_transfer_fraction(np.pi, 1.0, ratio), color="k", lw=1.2)
    for a, label in ((1, "H"), (2, "D"), (12, "C"), (207, "Pb")):
        value = float(energy_transfer_fraction(np.pi, 1.0, float(a)))
        ax_energy.plot(a, value, "o", color="tab:purple", ms=5)
        ax_energy.annotate(
            f"{label} {value:.3f}",
            (a, value),
            textcoords="offset points",
            xytext=(5, 3),
            fontsize=8,
        )
    ax_energy.set_xscale("log")
    ax_energy.set_xlim(0.01, 1000)
    ax_energy.set_ylim(0, 1.15)
    ax_energy.set_xlabel(r"質量比 $m_2/m_1$")
    ax_energy.set_ylabel("渡るエネルギーの割合")
    ax_energy.set_title(r"正面衝突（$\chi=180^\circ$）", fontsize=9)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_collision_video(
    batch: CollisionBatch, path: Path, fps: int = 30, body_radius: float = 0.18
) -> Path:
    """1 回の衝突を、実験室系・重心系・運動量の平面で同時に見る（実時間）。"""
    setup_matplotlib()
    q = batch.q[:, 0]
    v = batch.v[:, 0]
    cm = _center_of_mass(batch)[:, 0]
    dt = float(batch.t[1] - batch.t[0])
    frames = np.round(np.arange(0.0, batch.t[-1], 1.0 / fps) / dt).astype(int)

    fig, (ax_lab, ax_cm, ax_p) = plt.subplots(1, 3, figsize=(12.0, 4.4))
    views = []
    for ax, shift, title in (
        (ax_lab, np.zeros_like(cm), "実験室系"),
        (ax_cm, cm, "重心系"),
    ):
        ax.set_aspect("equal")
        ax.set_xlim(-6, 6)
        ax.set_ylim(-6, 6)
        ax.set_title(title)
        trails = [ax.plot([], [], color=c, lw=1)[0] for c in ("tab:blue", "tab:red")]
        bodies = [Circle((0, 0), body_radius, color=c) for c in ("tab:blue", "tab:red")]
        for body in bodies:
            ax.add_patch(body)
        views.append((shift, trails, bodies))
    (center,) = ax_lab.plot([], [], "k+", ms=9)

    total = batch.m1 * batch.speed
    circle_x, radius = momentum_circle(batch.m1, batch.m2, batch.speed)
    phi = np.linspace(0.0, 2.0 * np.pi, 300)
    # 衝突後に乗る円: 標的は実線、入射粒子（m1 v − p2'）は点線
    for middle, style in ((circle_x, "-"), (total - circle_x, ":")):
        ax_p.plot(
            (middle + radius * np.cos(phi)) / total,
            radius * np.sin(phi) / total,
            style,
            color="0.6",
        )
    ax_p.set_aspect("equal")
    ax_p.set_xlim(-0.45, 1.45)
    ax_p.set_ylim(-0.95, 0.95)
    ax_p.set_title("運動量（灰色: 衝突後に乗る円）")
    ax_p.set_xlabel(r"$p_x/m_1v$")
    ax_p.set_ylabel(r"$p_y/m_1v$")
    p1 = batch.m1 * v[:, 0:2] / total
    p2 = batch.m2 * v[:, 2:4] / total
    (path_line,) = ax_p.plot([], [], color="tab:red", lw=1)
    (dot1,) = ax_p.plot([], [], "o", color="tab:blue", ms=8, label="入射")
    (dot2,) = ax_p.plot([], [], "o", color="tab:red", ms=8, label="標的")
    ax_p.legend(loc="upper left", fontsize=9)
    clock = ax_lab.text(0.03, 0.94, "", transform=ax_lab.transAxes)
    fig.tight_layout()

    def update(i: int):
        for shift, trails, bodies in views:
            for j, (trail, body) in enumerate(zip(trails, bodies, strict=True)):
                position = q[:, 2 * j : 2 * j + 2] - shift
                trail.set_data(position[: i + 1, 0], position[: i + 1, 1])
                body.center = tuple(position[i])
        center.set_data([cm[i, 0]], [cm[i, 1]])
        path_line.set_data(p2[: i + 1, 0], p2[: i + 1, 1])
        dot1.set_data([p1[i, 0]], [p1[i, 1]])
        dot2.set_data([p2[i, 0]], [p2[i, 1]])
        clock.set_text(f"t = {batch.t[i]:.2f} s")
        return ()

    anim = animation.FuncAnimation(fig, update, frames=frames, blit=False)
    saved = save_animation(anim, path, fps=fps)
    plt.close(fig)
    return saved
