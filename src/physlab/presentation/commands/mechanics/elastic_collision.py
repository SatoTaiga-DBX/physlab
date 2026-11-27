"""`physlab run elastic-collision`: 力の形が違っても、衝突後の運動量は同じ円に乗る。

設計書: docs/design/007-collisions.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from physlab.application.mechanics.elastic_collision import (
    DEMO_POTENTIALS,
    circle_experiment,
    simulate_collisions,
)
from physlab.domain.mechanics.collision import max_deflection, momentum_circle
from physlab.infrastructure.rendering.mechanics.elastic_collision import (
    NAMES,
    render_collision_video,
    render_mass_ratio_figure,
    render_momentum_figure,
    render_trajectory_figure,
)
from physlab.infrastructure.storage import figure_path, output_path

NAME = "elastic-collision"
TITLE = "弾性衝突：力の形が違っても、衝突後の運動量は同じ円の上"
TOPIC = "mechanics"
RATIOS = (0.5, 1.0, 2.0)  # m2 / m1
IMPACTS = np.linspace(0.0, 1.2, 13)


def figures() -> list[Path]:
    """教科書の図 3 枚。"""
    cases = [(1.0, ratio, circle_experiment(1.0, ratio, IMPACTS)) for ratio in RATIOS]
    batches = {
        name: simulate_collisions(1.0, 2.0, potential, [0.6])
        for name, potential in DEMO_POTENTIALS.items()
    }
    points = [(ratio, results["gaussian"]) for _, ratio, results in cases]
    return [
        render_momentum_figure(cases, figure_path("col-momentum")),
        render_trajectory_figure(batches, figure_path("col-trajectories")),
        render_mass_ratio_figure(points, figure_path("col-mass-ratio")),
    ]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.add_argument("--m2", type=float, default=2.0, help="標的の質量（入射粒子は 1 kg）")
    p.add_argument("--impact", type=float, default=0.6, help="衝突パラメータ [m]")
    p.add_argument("--fps", type=int, default=30, help="動画のフレームレート")
    p.add_argument("--no-video", action="store_true", help="動画を作らず、図だけ作る")
    args = p.parse_args(argv)
    if args.m2 <= 0:
        p.error("--m2 は正の値")

    print(f"[{NAME}] m1 = 1 kg、m2 = {args.m2:g} kg、v = 1 m/s、ρ = {args.impact:g} m")
    print(f"  入射粒子の角の上限（保存則）: {np.degrees(max_deflection(1.0, args.m2)):.2f}°")
    center, radius = momentum_circle(1.0, args.m2, 1.0)
    results = circle_experiment(1.0, args.m2, [args.impact])
    print("  力の形ごとの結果（χ: 重心系、θ1・θ2: 実験室系）")
    for name, result in results.items():
        distance = np.hypot(result.p2[0, 0] - center, result.p2[0, 1])
        label = NAMES[name].replace("$", "")
        print(
            f"    {label:12s} χ={np.degrees(result.chi[0]):7.2f}°  "
            f"θ1={np.degrees(result.theta1[0]):7.2f}°  "
            f"θ2={np.degrees(result.theta2[0]):7.2f}°  "
            f"E2'/E1={result.transfer[0]:.4f}  円からのずれ {abs(distance - radius):.1e}"
        )

    print("図を作成中...")
    for path in figures():
        print(f"  -> {path}")
    if not args.no_video:
        print("動画を作成中...")
        batch = simulate_collisions(
            1.0, args.m2, DEMO_POTENTIALS["gaussian"], [args.impact], start=6.0, duration=12.0
        )
        video = render_collision_video(
            batch, output_path(TOPIC, "elastic_collision", ".mp4"), fps=args.fps
        )
        print(f"  -> {video}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
