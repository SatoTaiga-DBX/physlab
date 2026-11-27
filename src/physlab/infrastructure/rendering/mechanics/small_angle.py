"""微小角散乱の図。結果データを受け取って描くだけで、計算はしない。

設計書: docs/design/007-collisions.md
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import NullFormatter

from physlab.application.mechanics.small_angle import (
    CrossSectionTail,
    ErrorSweep,
    ImpulseDemo,
)
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_figure

NAMES = {
    "coulomb": r"$\alpha/r$",
    "inverse-square": r"$\alpha/r^2$",
    "inverse-fourth": r"$\alpha/r^4$",
    "gaussian": "ガウス型",
}
COLORS = ("tab:blue", "tab:orange", "tab:green", "tab:red")


def render_impulse_figure(demo: ImpulseDemo, path: Path) -> Path:
    """ほとんど直進する軌道と、横向きの力・運動量。直線とみなした近似と比べる。"""
    setup_matplotlib()
    fig, (ax_orbit, ax_force, ax_momentum) = plt.subplots(
        3, 1, figsize=(6.2, 5.4), constrained_layout=True, sharex=True
    )
    ax_orbit.plot(demo.x, demo.y, color="tab:blue", lw=1.6, label="実際の軌道")
    ax_orbit.axhline(demo.rho, color="k", ls="--", lw=1, label=r"直線 $y=\rho$")
    ax_orbit.set_ylim(0.4, 1.35)
    ax_orbit.set_ylabel("$y/a$")
    ax_orbit.set_title(
        f"偏向角: 数値解 {np.degrees(demo.exact_angle):.3f} 度、"
        f"力積の近似 {np.degrees(demo.approx_angle):.3f} 度",
        fontsize=10,
    )
    ax_orbit.legend(fontsize=8, loc="upper left")

    ax_force.plot(demo.x, demo.force_y, color="tab:blue", lw=1.6, label="実際の軌道の上")
    ax_force.plot(demo.straight_x, demo.straight_force_y, "k--", lw=1, label="直線の上")
    ax_force.set_ylabel("横向きの力 $F_y$")
    ax_force.legend(fontsize=8, loc="upper left")

    ax_momentum.plot(demo.x, demo.momentum_y, color="tab:blue", lw=1.6, label="実際の $p_y$")
    ax_momentum.axhline(demo.approx_momentum_y, color="k", ls="--", lw=1, label="力積の近似の終値")
    ax_momentum.set_ylabel("横向きの運動量 $p_y$")
    ax_momentum.set_xlabel("$x/a$")
    ax_momentum.legend(fontsize=8, loc="upper left")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_error_figure(
    sweeps: Mapping[str, ErrorSweep], tails: Sequence[CrossSectionTail], path: Path
) -> Path:
    """左: 力積の近似の相対誤差。右: 小さい角の断面積。どちらも両対数。"""
    setup_matplotlib()
    fig, (ax_error, ax_tail) = plt.subplots(1, 2, figsize=(6.5, 3.3), constrained_layout=True)
    for (name, sweep), color in zip(sweeps.items(), COLORS, strict=True):
        ax_error.plot(
            sweep.exact,
            sweep.relative_error,
            "o-",
            color=color,
            ms=3,
            lw=1,
            label=f"{NAMES[name]}（傾き {sweep.slope:.2f}）",
        )
    theta = np.array([4e-3, 0.5])
    ax_error.plot(theta, 0.9 * theta, color="0.6", lw=0.8, ls=":")
    ax_error.plot(theta, theta**2 / 12.0 * 0.5, color="0.6", lw=0.8, ls="--")
    ax_error.text(0.2, 0.55, "傾き 1", color="0.4", fontsize=7.5)
    ax_error.text(0.2, 6e-4, "傾き 2", color="0.4", fontsize=7.5)
    ax_error.set_xscale("log")
    ax_error.set_yscale("log")
    ax_error.set_xlabel(r"偏向角 $\theta$ [rad]")
    ax_error.set_ylabel("近似の相対誤差")
    ax_error.set_title("角が小さいほど正確", fontsize=10)
    ax_error.set_ylim(top=60)
    ax_error.legend(fontsize=7.5, loc="upper left")

    for tail, color in zip(tails, COLORS, strict=False):
        ax_tail.plot(
            tail.theta,
            tail.cross_section,
            color=color,
            lw=1.4,
            label=rf"$n={tail.n:g}$：傾き {tail.slope:.2f}（理論 {-2 - 2 / tail.n:g}）",
        )
    ax_tail.set_xscale("log")
    ax_tail.set_yscale("log")
    ax_tail.set_xticks(
        [2e-3, 5e-3, 1e-2, 2e-2, 5e-2], labels=["0.002", "0.005", "0.01", "0.02", "0.05"]
    )
    ax_tail.xaxis.set_minor_formatter(NullFormatter())
    ax_tail.set_xlabel(r"偏向角 $\theta$ [rad]")
    ax_tail.set_ylabel(r"$d\sigma/d\Omega$")
    ax_tail.set_title(r"$U=\alpha/r^n$ の小さい角の断面積", fontsize=10)
    ax_tail.set_ylim(top=ax_tail.get_ylim()[1] * 100)
    ax_tail.legend(fontsize=7.5, loc="upper right")
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved
