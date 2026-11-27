"""散乱の図と動画。結果データを受け取って描くだけで、計算はしない。

設計書: docs/design/007-collisions.md
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from physlab.application.mechanics.deflection import Orbits
from physlab.application.mechanics.scattering import (
    CrossSection,
    DeflectionCurve,
    LabCrossSection,
    branch_cross_section,
)
from physlab.domain.mechanics.decay import angle_density, max_lab_angle
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_animation, save_figure

RHO_MAP = plt.get_cmap("viridis")


@dataclass(frozen=True)
class CurveSeries:
    """偏向角の曲線と、時間積分で測った点。"""

    label: str
    curve: DeflectionCurve
    color: str
    orbit_rho: np.ndarray | None = None
    orbit_chi: np.ndarray | None = None


def _background(ax, width: float, extent: float) -> None:
    x = np.linspace(-extent, extent, 200)
    gx, gy = np.meshgrid(x, x)
    ax.imshow(
        np.exp(-(gx**2 + gy**2) / width**2),
        extent=(-extent, extent, -extent, extent),
        cmap="Oranges",
        vmin=0,
        vmax=1.6,
        origin="lower",
    )


def render_deflection_figure(
    orbits: Orbits,
    series: Sequence[CurveSeries],
    path: Path,
    energy_ratio: float,
    width: float = 1.0,
) -> Path:
    """左: 丘を通る粒子の束。右: 衝突パラメータと散乱角。"""
    setup_matplotlib()
    fig, (ax_orbit, ax_curve) = plt.subplots(1, 2, figsize=(6.5, 3.2), constrained_layout=True)
    extent = 4.5
    _background(ax_orbit, width, extent)
    norm = Normalize(0.0, float(np.max(np.abs(orbits.rho))))
    for j, rho in enumerate(orbits.rho):
        ax_orbit.plot(orbits.x[:, j], orbits.y[:, j], color=RHO_MAP(norm(abs(rho))), lw=0.8)
    ax_orbit.set_xlim(-extent, extent)
    ax_orbit.set_ylim(-extent, extent)
    ax_orbit.set_aspect("equal")
    ax_orbit.grid(False)
    ax_orbit.set_title(f"ガウス型の丘（$E/U_0={energy_ratio:g}$）", fontsize=10)
    ax_orbit.set_xlabel("$x/a$")
    ax_orbit.set_ylabel("$y/a$")
    fig.colorbar(
        ScalarMappable(norm=norm, cmap=RHO_MAP), ax=ax_orbit, label=r"$\rho/a$", shrink=0.8
    )

    for item in series:
        ax_curve.plot(
            item.curve.rho,
            np.degrees(item.curve.chi),
            color=item.color,
            lw=1.3,
            label=item.label,
        )
        if item.orbit_rho is not None:
            ax_curve.plot(item.orbit_rho, np.degrees(item.orbit_chi), "o", color=item.color, ms=3.5)
    ax_curve.set_xlim(0, 2.5)
    ax_curve.set_ylim(0, 245)
    ax_curve.set_yticks(range(0, 181, 45))
    ax_curve.set_xlabel(r"衝突パラメータ $\rho/a$")
    ax_curve.set_ylabel(r"散乱角 $\chi$ [度]")
    ax_curve.set_title("線: 求積、点: 運動方程式", fontsize=10)
    ax_curve.legend(fontsize=7.5, loc="upper right", ncol=2)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_cross_section_figure(
    sphere: CrossSection,
    hill: CrossSection,
    hill_curve: DeflectionCurve,
    hill_theory: np.ndarray,
    rainbow_angle: float,
    path: Path,
) -> Path:
    """乱数で打ち込んで数えた ``dσ/dΩ``。左: 剛体球、右: ガウス型の丘。"""
    setup_matplotlib()
    fig, (ax_sphere, ax_hill) = plt.subplots(1, 2, figsize=(6.5, 3.1), constrained_layout=True)
    ax_sphere.errorbar(
        np.degrees(sphere.centers), sphere.measured, yerr=sphere.error, fmt="o", ms=3.5
    )
    ax_sphere.axhline(0.25, color="k", lw=1, label=r"理論 $a^2/4$")
    ax_sphere.set_ylim(0, 0.4)
    ax_sphere.set_title("剛体球：どの向きにも同じ", fontsize=10)
    ax_sphere.set_ylabel(r"$d\sigma/d\Omega\ /\,a^2$")
    ax_sphere.legend(fontsize=8)

    used = hill.counts > 0
    ax_hill.errorbar(
        np.degrees(hill.centers[used]),
        hill.measured[used],
        yerr=hill.error[used],
        fmt="o",
        ms=3.5,
        color="tab:blue",
        label="数えた値",
    )
    ax_hill.step(
        np.degrees(hill.edges),
        np.append(hill_theory, hill_theory[-1]),
        where="post",
        color="0.6",
        lw=1,
        label="理論（区間の平均）",
    )
    chi = np.linspace(np.radians(2.0), rainbow_angle - np.radians(0.05), 300)
    ax_hill.plot(
        np.degrees(chi),
        branch_cross_section(hill_curve, chi),
        color="k",
        lw=1,
        label="理論（点ごと）",
    )
    ax_hill.axvline(np.degrees(rainbow_angle), color="tab:red", ls="--", lw=1)
    ax_hill.text(
        np.degrees(rainbow_angle) - 0.8,
        0.5,
        "虹の角",
        color="tab:red",
        fontsize=8,
        ha="right",
    )
    ax_hill.set_yscale("log")
    ax_hill.set_ylim(0.2, 200)
    ax_hill.set_xlim(0, np.degrees(hill.edges[-1]))
    ax_hill.set_title("ガウス型の丘：虹の角に集まる", fontsize=10)
    ax_hill.legend(fontsize=7.5, loc="upper right")
    for ax in (ax_sphere, ax_hill):
        ax.set_xlabel(r"散乱角 $\chi$ [度]")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_lab_figure(results: Sequence[LabCrossSection], path: Path) -> Path:
    """剛体球どうしの衝突を実験室系で見た ``dσ/dΩ1``。"""
    setup_matplotlib()
    fig, axes = plt.subplots(1, len(results), figsize=(6.5, 2.6), constrained_layout=True)
    for ax, result in zip(axes, results, strict=True):
        counted = result.counted
        k = result.mass_ratio
        used = counted.counts > 0
        ax.errorbar(
            np.degrees(counted.centers[used]),
            counted.measured[used],
            yerr=counted.error[used],
            fmt="o",
            ms=3,
            label="数えた値",
        )
        theta = np.linspace(1e-3, max_lab_angle(k) * (1 - 1e-4), 600)
        if k == 1.0:
            theta = theta[theta < np.pi / 2]
        theory = np.pi * angle_density(theta, k) / (2.0 * np.pi * np.sin(theta))
        ax.plot(np.degrees(theta), theory, color="k", lw=1, label="理論")
        ax.set_title(f"$m_1/m_2={k:g}$", fontsize=10)
        ax.set_xlabel(r"$\theta_1$ [度]")
        ax.set_xlim(0, np.degrees(max_lab_angle(k)))
        if k > 1.0:  # 上限の角で無限大になるので対数軸
            ax.set_yscale("log")
            ax.set_ylim(5, 300)
        else:
            ax.set_ylim(0, 1.3 * float(np.max(theory)))
    axes[0].set_ylabel(r"$d\sigma/d\Omega_1\ /\,a^2$")
    axes[0].legend(fontsize=7.5)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_beam_video(
    orbits: Orbits,
    chi: np.ndarray,
    curve: DeflectionCurve,
    path: Path,
    width: float = 1.0,
    fps: int = 30,
) -> Path:
    """左: 粒子の束が丘を通る。右: 出ていった粒子の ``(ρ, χ)`` が曲線の上に現れる（実時間）。"""
    setup_matplotlib()
    t = orbits.t
    dt = float(t[1] - t[0])
    frames = np.round(np.arange(0.0, t[-1], 1.0 / fps) / dt).astype(int)
    extent = 6.0

    fig, (ax_orbit, ax_curve) = plt.subplots(1, 2, figsize=(10.0, 5.0))
    _background(ax_orbit, width, extent)
    ax_orbit.set_xlim(-extent, extent)
    ax_orbit.set_ylim(-extent, extent)
    ax_orbit.set_aspect("equal")
    ax_orbit.grid(False)
    ax_orbit.set_title("同じ速さの粒子の束")
    norm = Normalize(0.0, float(np.max(np.abs(orbits.rho))))
    colors = RHO_MAP(norm(np.abs(orbits.rho)))
    trails = [ax_orbit.plot([], [], color=c, lw=0.8)[0] for c in colors]
    heads = ax_orbit.scatter(orbits.x[0], orbits.y[0], c=colors, s=14)

    ax_curve.plot(curve.rho, np.degrees(curve.chi), color="0.7", lw=1.2)
    limit = float(np.max(np.abs(orbits.rho)))
    ax_curve.set_xlim(0, limit * 1.02)
    ax_curve.set_ylim(0, np.degrees(curve.chi.max()) * 1.15)
    ax_curve.set_xlabel(r"$|\rho|/a$")
    ax_curve.set_ylabel(r"散乱角 $\chi$ [度]")
    ax_curve.set_title(r"力の中心を通り過ぎた粒子の $(\rho,\ \chi)$")
    found = ax_curve.scatter([], [], s=22)
    clock = ax_orbit.text(0.02, 0.95, "", transform=ax_orbit.transAxes)
    fig.tight_layout()
    radius = np.hypot(orbits.x, orbits.y)
    outgoing = orbits.x * orbits.vx + orbits.y * orbits.vy > 0

    def update(i: int):
        for j, trail in enumerate(trails):
            trail.set_data(orbits.x[: i + 1, j], orbits.y[: i + 1, j])
        heads.set_offsets(np.column_stack([orbits.x[i], orbits.y[i]]))
        done = outgoing[i] & (radius[i] > 4.0)
        found.set_offsets(
            np.column_stack([np.abs(orbits.rho[done]), np.degrees(chi[done])]).reshape(-1, 2)
        )
        found.set_color(colors[done])
        clock.set_text(f"t = {t[i]:.2f} s")
        return ()

    anim = animation.FuncAnimation(fig, update, frames=frames, blit=False)
    saved = save_animation(anim, path, fps=fps)
    plt.close(fig)
    return saved
