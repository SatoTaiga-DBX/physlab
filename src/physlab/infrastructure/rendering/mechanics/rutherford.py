"""ラザフォード散乱の図と動画。結果データを受け取って描くだけで、計算はしない。

設計書: docs/design/007-collisions.md
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
from matplotlib.patches import Circle

from physlab.application.mechanics.deflection import Orbits
from physlab.application.mechanics.rutherford import EnergyTransfer
from physlab.application.mechanics.scattering import CrossSection, DeflectionCurve
from physlab.domain.mechanics.scattering import rutherford_angle, rutherford_cross_section
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_animation, save_figure

BALL_COLORS = ("tab:blue", "tab:orange", "tab:green")


def _frame(ax, head_on: float, extent: float) -> None:
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_aspect("equal")
    ax.add_patch(Circle((0, 0), head_on, fill=False, ls="--", color="0.5", lw=0.9))
    ax.plot([0], [0], "o", color="tab:red", ms=6)


def render_orbit_figure(repulsive: Orbits, attractive: Orbits, head_on: float, path: Path) -> Path:
    """左: 斥力の軌道の束。右: 同じ ρ の斥力と引力。曲がる角の大きさは同じ。"""
    setup_matplotlib()
    fig, (ax_beam, ax_pair) = plt.subplots(1, 2, figsize=(6.5, 3.3), constrained_layout=True)
    _frame(ax_beam, head_on, 5.0)
    for j in range(repulsive.rho.size):
        ax_beam.plot(repulsive.x[:, j], repulsive.y[:, j], color="tab:blue", lw=0.8)
    ax_beam.set_title("斥力：中心の裏側には影ができる", fontsize=9)
    ax_beam.set_ylabel("$y/d$")

    _frame(ax_pair, head_on, 5.0)
    pick = attractive.rho.size // 2
    rho = attractive.rho[pick]
    j = int(np.argmin(np.abs(repulsive.rho - rho)))
    ax_pair.plot(repulsive.x[:, j], repulsive.y[:, j], color="tab:blue", lw=1.5, label="斥力")
    ax_pair.plot(
        attractive.x[:, pick], attractive.y[:, pick], color="tab:red", lw=1.5, label="引力"
    )
    ax_pair.axhline(rho, color="k", lw=0.7, ls=":")
    ax_pair.text(-4.8, rho + 0.2, rf"$\rho={rho:g}\,d$", fontsize=8)
    chi = np.degrees(rutherford_angle(rho, head_on))
    ax_pair.set_title(f"曲がる角はどちらも {chi:.1f}°", fontsize=9)
    ax_pair.legend(fontsize=8, loc="lower right")
    for ax in (ax_beam, ax_pair):
        ax.set_xlabel("$x/d$")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_cross_section_figure(
    counted: CrossSection, transfer: EnergyTransfer, head_on: float, path: Path
) -> Path:
    """左: 散乱角の分布。右: 標的に渡るエネルギーの分布。点は数えた値、線は公式。"""
    setup_matplotlib()
    fig, (ax_angle, ax_energy) = plt.subplots(1, 2, figsize=(6.5, 3.2), constrained_layout=True)
    ax_angle.errorbar(
        np.degrees(counted.centers),
        counted.measured,
        yerr=counted.error,
        fmt="o",
        ms=3.5,
        label="数えた値",
    )
    chi = np.linspace(counted.edges[0], counted.edges[-1], 300)
    ax_angle.plot(
        np.degrees(chi),
        rutherford_cross_section(chi, head_on),
        color="k",
        lw=1,
        label=r"$(d/4)^2/\sin^4(\chi/2)$",
    )
    ax_angle.set_yscale("log")
    ax_angle.set_xlabel(r"散乱角 $\chi$ [度]")
    ax_angle.set_ylabel(r"$d\sigma/d\Omega\ /\,d^2$")
    ax_angle.set_title("向きの分布", fontsize=10)
    ax_angle.legend(fontsize=8)

    ratio = transfer.centers / transfer.maximum
    ax_energy.errorbar(
        ratio,
        transfer.measured,
        yerr=transfer.error,
        fmt="o",
        ms=3.5,
        label="数えた値",
    )
    ax_energy.plot(ratio, transfer.expected, color="k", lw=1, label=r"公式 $\propto 1/\epsilon^2$")
    ax_energy.set_xscale("log")
    ax_energy.set_yscale("log")
    ax_energy.set_xlabel(r"渡すエネルギー $\epsilon/\epsilon_{\max}$")
    ax_energy.set_ylabel(r"$d\sigma/d\epsilon$")
    ax_energy.set_title(
        rf"標的に渡るエネルギー（$m_2/m_1={transfer.m2 / transfer.m1:g}$）", fontsize=10
    )
    ax_energy.legend(fontsize=8)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_charge_size_figure(
    curves: Mapping[float, DeflectionCurve],
    counted: Mapping[float, CrossSection],
    head_on: float,
    path: Path,
) -> Path:
    """電荷を半径 R の球に広げる。左: 散乱角、右: 数えた断面積。破線は点電荷。"""
    setup_matplotlib()
    fig, (ax_curve, ax_count) = plt.subplots(1, 2, figsize=(6.5, 3.2), constrained_layout=True)
    rho = next(iter(curves.values())).rho
    ax_curve.plot(
        rho, np.degrees(rutherford_angle(rho, head_on)), "k--", lw=1.2, label="点電荷", zorder=5
    )
    for (ratio, curve), color in zip(curves.items(), BALL_COLORS, strict=True):
        ax_curve.plot(
            curve.rho, np.degrees(curve.chi), color=color, lw=2.0, label=f"$R={ratio:g}\\,d$"
        )
    ax_curve.set_xlim(0, 4)
    ax_curve.set_ylim(0, 185)
    ax_curve.set_yticks(range(0, 181, 45))
    ax_curve.set_xlabel(r"衝突パラメータ $\rho/d$")
    ax_curve.set_ylabel(r"散乱角 $\chi$ [度]")
    ax_curve.legend(fontsize=7.5, loc="upper right")

    chi = np.linspace(np.radians(20.0), np.pi, 300)
    ax_count.plot(
        np.degrees(chi), rutherford_cross_section(chi, head_on), "k--", lw=1.2, label="点電荷"
    )
    for (ratio, result), color, marker in zip(
        counted.items(), BALL_COLORS, ("o", "s", "^"), strict=True
    ):
        used = result.counts > 0
        ax_count.errorbar(
            np.degrees(result.centers[used]),
            result.measured[used],
            yerr=result.error[used],
            fmt=marker,
            ms=3.5,
            color=color,
            label=f"$R={ratio:g}\\,d$",
        )
    ax_count.set_yscale("log")
    ax_count.set_xlim(0, 180)
    ax_count.set_xticks(range(0, 181, 45))
    ax_count.set_xlabel(r"散乱角 $\chi$ [度]")
    ax_count.set_ylabel(r"$d\sigma/d\Omega\ /\,d^2$")
    ax_count.legend(fontsize=7.5, loc="upper right")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_rutherford_video(
    orbits: Orbits, head_on: float, path: Path, fps: int = 30, extent: float = 6.0
) -> Path:
    """衝突パラメータの違う粒子が、原子核に近づいて曲がる（実時間）。"""
    setup_matplotlib()
    t = orbits.t
    dt = float(t[1] - t[0])
    frames = np.round(np.arange(0.0, t[-1], 1.0 / fps) / dt).astype(int)
    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    _frame(ax, head_on, extent)
    ax.set_title("斥力の中心に打ち込む（破線: 半径 d の円）")
    ax.set_xlabel("$x/d$")
    ax.set_ylabel("$y/d$")
    cmap = plt.get_cmap("viridis")
    scale = float(np.max(np.abs(orbits.rho)))
    colors = [cmap(abs(r) / scale) for r in orbits.rho]
    trails = [ax.plot([], [], color=c, lw=0.9)[0] for c in colors]
    heads = ax.scatter(orbits.x[0], orbits.y[0], c=colors, s=18)
    clock = ax.text(0.02, 0.96, "", transform=ax.transAxes)
    fig.tight_layout()

    def update(i: int):
        for j, trail in enumerate(trails):
            trail.set_data(orbits.x[: i + 1, j], orbits.y[: i + 1, j])
        heads.set_offsets(np.column_stack([orbits.x[i], orbits.y[i]]))
        clock.set_text(f"t = {t[i]:.2f} s")
        return ()

    anim = animation.FuncAnimation(fig, update, frames=frames, blit=False)
    saved = save_animation(anim, path, fps=fps)
    plt.close(fig)
    return saved
