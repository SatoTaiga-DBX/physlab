"""力学的相似の図。結果データを受け取って描くだけで、計算はしない。

設計書: docs/design/006-conservation-laws.md
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import NullFormatter

from physlab.application.mechanics.similarity import Collapse, KeplerSweep, PeriodSweep, Virial
from physlab.infrastructure.rendering.style import setup_matplotlib
from physlab.infrastructure.rendering.writers import save_figure

K_COLORS = {1.0: "tab:purple", 2.0: "tab:blue", 4.0: "tab:green", 6.0: "tab:red"}


def render_period_figure(sweep: PeriodSweep, kepler: KeplerSweep, path: Path) -> Path:
    """周期と大きさの両対数グラフ。傾きが 1 − k/2（万有引力は 3/2）。"""
    setup_matplotlib()
    fig, ax = plt.subplots(figsize=(6.2, 3.8), constrained_layout=True)
    fine = np.geomspace(sweep.amplitudes[0], sweep.amplitudes[-1], 100)
    for k, measured in sweep.measured.items():
        coefficient = sweep.exact[k][0] / sweep.amplitudes[0] ** (1 - k / 2)
        ax.plot(fine, coefficient * fine ** (1 - k / 2), color=K_COLORS[k], lw=1)
        ax.plot(
            sweep.amplitudes,
            measured,
            "o",
            color=K_COLORS[k],
            ms=5,
            label=f"$k={k:g}$：傾き {sweep.slopes[k]:+.3f}（理論 {1 - k / 2:+g}）",
        )
    ax.plot(
        kepler.semi_major,
        kepler.measured,
        "s",
        color="k",
        ms=5,
        label=f"万有引力 $k=-1$：傾き {kepler.slope:+.3f}（理論 +1.5）",
    )
    a = np.geomspace(kepler.semi_major.min(), kepler.semi_major.max(), 50)
    ax.plot(a, kepler.exact[0] * (a / kepler.semi_major[0]) ** 1.5, color="k", lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("大きさ（振幅 $A$、万有引力は長半径 $a$）")
    ax.set_ylabel("周期 $T$")
    ax.set_title("周期は大きさの $1-k/2$ 乗に比例する（点: 数値解、線: 理論）", fontsize=10)
    ticks = [0.25, 0.5, 1.0, 2.0, 4.0]  # 対数軸の目盛りは詰まらないように固定する
    ax.set_xticks(ticks, labels=[f"{v:g}" for v in ticks])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(top=ax.get_ylim()[1] * 15)  # 凡例の余白
    ax.legend(fontsize=7.5, loc="upper center", ncol=2)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_collapse_figure(result: Collapse, path: Path) -> Path:
    """左: 振幅の違う運動は周期も違う。右: 大きさと時間を相似に直すと 1 本に重なる。"""
    setup_matplotlib()
    fig, (ax_raw, ax_scaled) = plt.subplots(1, 2, figsize=(6.4, 2.8), constrained_layout=True)
    widths = (5.0, 2.8, 1.2)
    colors = ("0.7", "tab:orange", "k")
    for a, t, x, st, sx, width, color in zip(
        result.amplitudes,
        result.t,
        result.x,
        result.scaled_t,
        result.scaled_x,
        widths,
        colors,
        strict=True,
    ):
        ax_raw.plot(t, x, color=color, lw=min(width, 1.8), label=f"$A={a:g}$")
        ax_scaled.plot(st, sx, color=color, lw=width, label=f"$A={a:g}$")
    ax_raw.set_xlabel("時間 $t$")
    ax_raw.set_ylabel("$x$")
    ax_raw.set_title(f"そのまま（$U=|x|^{{{result.k:g}}}/{result.k:g}$）", fontsize=10)
    ax_scaled.set_xlabel(r"$t\,/\,A^{1-k/2}$")
    ax_scaled.set_ylabel("$x/A$")
    ax_scaled.set_title("相似に直すと重なる", fontsize=10)
    ax_scaled.set_xlim(0, result.scaled_t[0][-1])
    for ax in (ax_raw, ax_scaled):
        low, high = ax.get_ylim()
        ax.set_ylim(low, high + 0.35 * (high - low))
        ax.legend(fontsize=7.5, loc="upper center", ncol=3)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved


def render_virial_figure(result: Virial, path: Path) -> Path:
    """時刻 0 から t までの平均の比 2⟨T⟩/⟨U⟩ は、k に近づく（ビリアル定理）。"""
    setup_matplotlib()
    fig, ax = plt.subplots(figsize=(6.2, 3.0), constrained_layout=True)
    for k in result.ks:
        ax.plot(result.periods[k], result.running[k], color=K_COLORS[k], lw=1.4, label=f"$k={k:g}$")
        ax.axhline(k, color=K_COLORS[k], ls="--", lw=0.9)
    ax.set_xlim(0, result.periods[result.ks[0]][-1])
    ax.set_ylim(0, 9)
    ax.set_xlabel("時間（周期の数 $t/T$）")
    ax.set_ylabel(r"$2\langle T\rangle/\langle U\rangle$")
    ax.set_title("時間平均の比は $k$ に落ち着く（破線）", fontsize=10)
    ax.legend(fontsize=8, loc="upper right", ncol=4)
    saved = save_figure(fig, path)
    plt.close(fig)
    return saved
