"""対称性と保存量の図と動画。結果データ（application）を受け取って描くだけで、計算はしない。

設計書: docs/design/006-conservation-laws.md
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt

from physlab.application.mechanics.symmetry import ClusterRun, IntegratorComparison, RodRun
from physlab.domain.mechanics.cluster import (
    as_vectors,
    center_of_mass,
    external_energy,
    kinetic_energy,
    spring_energy,
)
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_animation, save_figure

FIELD_STYLES = {  # 外場の名前, 色, 線幅（重なっても見分けられるように太さを変える）
    "none": ("外場なし", "0.6", 5.0),
    "gravity": ("一様重力（下向き）", "tab:blue", 2.4),
    "trap": ("原点へ引く外場", "tab:orange", 1.6),
    "pulsing_trap": ("強さが脈打つ外場", "tab:red", 1.0),
}
PARTICLE_COLORS = ("tab:blue", "tab:orange", "tab:green")


def render_symmetry_figure(runs: dict[str, ClusterRun], path: Path) -> Path:
    """E, Px, Py, L の時間変化。点は「対称性の破れ」から力を積分して予測した値。"""
    setup_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(6.4, 4.6), sharex=True, constrained_layout=True)
    panels = [
        (axes[0, 0], "エネルギー $E$", lambda r: (r.energy, r.predicted_energy)),
        (axes[0, 1], "運動量 $P_x$", lambda r: (r.momentum[:, 0], r.predicted_momentum[:, 0])),
        (axes[1, 0], "運動量 $P_y$", lambda r: (r.momentum[:, 1], r.predicted_momentum[:, 1])),
        (axes[1, 1], "角運動量 $L$", lambda r: (r.angular, r.predicted_angular)),
    ]
    for ax, title, pick in panels:
        for field, run in runs.items():
            label, color, width = FIELD_STYLES[field]
            value, predicted = pick(run)
            ax.plot(run.t, value, color=color, lw=width, label=label)
            if np.ptp(value) > 0.1:  # 変化する量にだけ予測の点を重ねる
                step = run.t.size // 12
                ax.plot(
                    run.t[step::step],
                    predicted[step::step],
                    "o",
                    ms=3.5,
                    mfc="white",
                    mec="k",
                    mew=0.8,
                    zorder=5,
                )
        ax.set_title(title, fontsize=10)
    axes[1, 0].set_ylim(-9, 7)
    axes[1, 0].text(
        0.97,
        0.04,
        "重力では直線的に減り続ける",
        transform=axes[1, 0].transAxes,
        ha="right",
        fontsize=7.5,
        color="tab:blue",
    )
    for ax in axes[1]:
        ax.set_xlabel("時間 $t$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=4, fontsize=8)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_integrator_figure(c: IntegratorComparison, path: Path) -> Path:
    """自由な系（全部保存するはず）で、速度ベルレと RK4 の保存量のずれを比べる。"""
    setup_matplotlib()
    fig, axes = plt.subplots(1, 3, figsize=(6.4, 2.5), sharey=True, constrained_layout=True)
    panels = [
        (axes[0], "エネルギー $|E-E_0|$", c.energy_error),
        (axes[1], "運動量 $|P-P_0|$", c.momentum_error),
        (axes[2], "角運動量 $|L-L_0|$", c.angular_error),
    ]
    for ax, title, errors in panels:
        for method, color, label in (
            ("verlet", "tab:green", "速度ベルレ"),
            ("rk4", "tab:orange", "RK4"),
        ):
            ax.plot(c.t, np.maximum(errors[method], 1e-17), color=color, lw=0.9, label=label)
        ax.set_yscale("log")
        ax.set_ylim(1e-17, 1e-2)
        ax.set_title(title, fontsize=9.5)
        ax.set_xlabel("時間 $t$")
    axes[0].set_ylabel("ずれ")
    axes[2].legend(fontsize=8, loc="center right")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_rod_figure(rod: RodRun, path: Path) -> Path:
    """回転する棒の上の質点: E = pṙ − L は一定、T は増える（増分はモーターの仕事）。"""
    setup_matplotlib()
    fig, ax = plt.subplots(figsize=(5.6, 2.8), constrained_layout=True)
    ax.plot(rod.t, rod.kinetic, color="tab:red", lw=2, label="運動エネルギー $T$")
    ax.plot(
        rod.t,
        rod.kinetic - rod.motor_work,
        color="0.3",
        lw=4,
        alpha=0.4,
        label="$T$ − モーターの仕事",
    )
    ax.plot(
        rod.t,
        rod.energy_function,
        color="tab:blue",
        lw=2,
        label=r"$E=p_r\dot r-\mathcal{L}$",
    )
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xlabel("時間 $t$")
    ax.set_ylabel("エネルギー")
    ax.set_title("$E$ は保存するが、$T$ は保存しない", fontsize=10)
    ax.legend(fontsize=8.5, loc="upper left")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_cluster_video(run: ClusterRun, path: Path, fps: int = 30) -> Path:
    """自由な 3 質点。個々の運動量やエネルギーの内訳は揺れ動くが、合計は一定。"""
    setup_matplotlib()
    p = run.params
    r = as_vectors(run.q)
    masses = np.asarray(p.masses)
    particle_px = masses * as_vectors(run.v)[:, :, 0]
    kinetic = kinetic_energy(run.v, p)
    potential = spring_energy(run.q, p) + external_energy(run.t, run.q, p)
    centre = center_of_mass(run.q, p)
    dt = float(run.t[1] - run.t[0])
    idx = np.round(np.arange(0.0, run.t[-1], 1.0 / fps) / dt).astype(int)

    fig = plt.figure(figsize=(10, 4.8))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.1, 1])
    ax_s = fig.add_subplot(grid[:, 0])
    ax_p = fig.add_subplot(grid[0, 1])
    ax_e = fig.add_subplot(grid[1, 1])

    pad = 0.4
    ax_s.set_xlim(r[..., 0].min() - pad, r[..., 0].max() + pad)
    ax_s.set_ylim(r[..., 1].min() - pad, r[..., 1].max() + pad)
    ax_s.set_aspect("equal")
    ax_s.set_title("ばねでつないだ 3 質点（外場なし）")
    (springs,) = ax_s.plot([], [], "-", color="0.6", lw=1)
    dots = [
        ax_s.plot([], [], "o", color=c, ms=8 + 5 * m)[0]
        for c, m in zip(PARTICLE_COLORS, masses, strict=True)
    ]
    (cm_trail,) = ax_s.plot([], [], ":", color="k", lw=1.2)
    (cm_dot,) = ax_s.plot([], [], "k+", ms=11, mew=2, label="重心")
    ax_s.legend(loc="upper right")
    time_text = ax_s.text(0.03, 0.04, "", transform=ax_s.transAxes)

    for a, color in enumerate(PARTICLE_COLORS):
        ax_p.plot(run.t, particle_px[:, a], color=color, lw=0.8, alpha=0.25)
    ax_p.plot(run.t, run.momentum[:, 0], color="k", lw=0.8, alpha=0.25)
    p_lines = [ax_p.plot([], [], color=c, lw=1.2)[0] for c in PARTICLE_COLORS]
    (p_total,) = ax_p.plot([], [], color="k", lw=2.5, label="合計 $P_x$（一定）")
    ax_p.set_xlim(run.t[0], run.t[-1])
    ax_p.set_title("$x$ 方向の運動量：色は各質点")
    ax_p.legend(loc="upper right", fontsize=8)

    ax_e.set_xlim(run.t[0], run.t[-1])
    ax_e.set_ylim(0, 1.25 * run.energy.max())
    (k_line,) = ax_e.plot([], [], color="tab:red", lw=1.2, label="運動エネルギー")
    (u_line,) = ax_e.plot([], [], color="tab:purple", lw=1.2, label="ばね")
    (e_line,) = ax_e.plot([], [], color="k", lw=2.5, label="合計 $E$（一定）")
    ax_e.set_xlabel("時間 $t$")
    ax_e.set_title("エネルギー")
    ax_e.legend(loc="upper right", fontsize=8, ncol=3)
    fig.tight_layout()

    loop = [0, 1, 2, 0]

    def update(i: int):
        springs.set_data(r[i, loop, 0], r[i, loop, 1])
        for a, dot in enumerate(dots):
            dot.set_data([r[i, a, 0]], [r[i, a, 1]])
        cm_trail.set_data(centre[: i + 1, 0], centre[: i + 1, 1])
        cm_dot.set_data([centre[i, 0]], [centre[i, 1]])
        for a, line in enumerate(p_lines):
            line.set_data(run.t[: i + 1], particle_px[: i + 1, a])
        p_total.set_data(run.t[: i + 1], run.momentum[: i + 1, 0])
        k_line.set_data(run.t[: i + 1], kinetic[: i + 1])
        u_line.set_data(run.t[: i + 1], potential[: i + 1])
        e_line.set_data(run.t[: i + 1], run.energy[: i + 1])
        time_text.set_text(f"t = {run.t[i]:5.2f}")
        return (
            springs,
            *dots,
            cm_trail,
            cm_dot,
            *p_lines,
            p_total,
            k_line,
            u_line,
            e_line,
            time_text,
        )

    anim = animation.FuncAnimation(fig, update, frames=idx, blit=True)
    saved = save_animation(anim, path, fps=fps)
    plt.close(fig)
    return saved
