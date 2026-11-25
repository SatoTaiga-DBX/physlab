"""`physlab run least-action`: 最小作用の原理（ハミルトンの原理）を、経路の族と作用の地形で見る。

設計書: docs/design/005-analytical-mechanics-ch1.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from physlab.application.mechanics.least_action import (
    action_landscape,
    relax_path,
    scan_action,
    scan_curvature_vs_duration,
)
from physlab.domain.mechanics.lagrangian import harmonic_oscillator, uniform_gravity
from physlab.infrastructure.rendering.mechanics.least_action import (
    render_curvature_figure,
    render_landscape_figure,
    render_path_family_figure,
    render_relaxation_figure,
    render_relaxation_video,
    render_second_variation_figure,
)
from physlab.infrastructure.storage import figure_path, output_path

NAME = "least-action"
TITLE = "最小作用の原理：経路の作用、最小か鞍点か、作用を下げると真の運動に収束する"
TOPIC = "mechanics"

GRAVITY = uniform_gravity(mass=1.0, g=9.81)
OSCILLATOR = harmonic_oscillator(mass=1.0, omega=2.0 * np.pi)  # π/ω = 0.5 s
TAU_SHORT, TAU_LONG = 0.3, 0.8  # 調和振動子で、π/ω より短い／長い時間
Q_START, Q_END = 1.0, 0.0  # 調和振動子の両端の位置
TAU_GRAVITY = 2.0  # 教科書の一様重力の例の時間の長さ


def figures() -> list[Path]:
    """教科書の図 5 枚（``outputs/figures/la-*.pdf``）。常に教科書の条件で作る。"""
    tau = TAU_GRAVITY
    eps = np.linspace(-2.0, 2.0, 21)
    eps_osc = np.linspace(-1.0, 1.0, 21)

    family = scan_action(GRAVITY, 0.0, 0.0, tau, 0.0, eps)
    short = scan_action(OSCILLATOR, 0.0, Q_START, TAU_SHORT, Q_END, eps_osc, mode=1)
    long1 = scan_action(OSCILLATOR, 0.0, Q_START, TAU_LONG, Q_END, eps_osc, mode=1)
    long2 = scan_action(OSCILLATOR, 0.0, Q_START, TAU_LONG, Q_END, eps_osc, mode=2)
    curvature = scan_curvature_vs_duration(OSCILLATOR, np.linspace(0.2, 1.6, 141))
    bowl = action_landscape(OSCILLATOR, 0.0, Q_START, TAU_SHORT, Q_END)
    saddle = action_landscape(OSCILLATOR, 0.0, Q_START, TAU_LONG, Q_END)
    relaxation = relax_path(GRAVITY, 0.0, 0.0, tau, 0.0)

    return [
        render_path_family_figure(family, figure_path("la-path-family")),
        render_second_variation_figure(short, long1, long2, figure_path("la-second-variation")),
        render_curvature_figure(curvature, figure_path("la-curvature")),
        render_landscape_figure(bowl, saddle, figure_path("la-landscape")),
        render_relaxation_figure(relaxation, figure_path("la-relaxation")),
    ]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.add_argument(
        "--duration", type=float, default=TAU_GRAVITY, help="一様重力の例の時間の長さ τ [s]"
    )
    p.add_argument("--fps", type=int, default=24, help="動画のフレームレート")
    p.add_argument("--no-video", action="store_true", help="動画を作らず、図だけ作る")
    args = p.parse_args(argv)

    tau = args.duration
    print(f"[{NAME}] 一様重力（m=1 kg, g=9.81）で (t, q) = (0, 0) → ({tau:g}, 0)")
    scan = scan_action(GRAVITY, 0.0, 0.0, tau, 0.0, np.linspace(-2.0, 2.0, 21))
    exact = -(9.81**2) * tau**3 / 24  # 解析解 −g²τ³/24
    print(f"  真の運動の作用 S* = {scan.action_true:.4f} J·s （解析解 {exact:.4f}）")
    print(
        f"  S(ε) = S* + c ε²: c = {scan.fit_curvature:.4f}（公式 {scan.theory_curvature:.4f}）、"
        f"1 次の係数 {scan.fit_linear:.1e}"
    )
    for tau_osc in (TAU_SHORT, TAU_LONG):
        cs = [OSCILLATOR.action_curvature(tau_osc, n) for n in (1, 2)]
        print(f"  調和振動子 τ={tau_osc} s: c1={cs[0]:+.3f}, c2={cs[1]:+.3f}（π/ω = 0.5 s）")

    print("図を作成中...")
    for path in figures():  # 引数を変えても、教科書の図は既定の条件のまま
        print(f"  -> {path}")

    if not args.no_video:
        print("動画を作成中...")
        relaxation = relax_path(GRAVITY, 0.0, 0.0, tau, 0.0)
        print(
            f"  作用 {relaxation.action[0]:.1f} → {relaxation.action[-1]:.3f}、"
            f"最大のずれ {relaxation.max_error[0]:.2f} → {relaxation.max_error[-1]:.1e} m"
        )
        video = render_relaxation_video(
            relaxation, output_path(TOPIC, "least_action_relaxation", ".mp4"), fps=args.fps
        )
        print(f"  -> {video}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
