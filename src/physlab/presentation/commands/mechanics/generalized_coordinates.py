"""`physlab run generalized-coordinates`: 一般座標と配位空間（単振り子と二重振り子）。

設計書: docs/design/005-analytical-mechanics-ch1.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from physlab.application.mechanics.generalized_coordinates import (
    measure_normal_mode_period,
    same_position_different_velocity,
    simulate_double_pendulum,
)
from physlab.domain.mechanics.double_pendulum import DoublePendulumParams
from physlab.infrastructure.rendering.mechanics.generalized_coordinates import (
    render_double_pendulum_figure,
    render_double_pendulum_video,
    render_same_coordinate_figure,
)
from physlab.infrastructure.storage import figure_path, output_path

NAME = "generalized-coordinates"
TITLE = "一般座標：位置だけでは運動は決まらない（単振り子）、1 点が系全体を表す（二重振り子）"
TOPIC = "mechanics"


def figures() -> list[Path]:
    """教科書の図 2 枚（``outputs/figures/gc-*.pdf``）。常に教科書の条件で作る。"""
    runs = same_position_different_velocity(np.deg2rad(60.0), np.array([0.0, 2.0, -3.0, 4.0]))
    double = simulate_double_pendulum(DoublePendulumParams(), 10.0)
    return [
        render_same_coordinate_figure(runs, figure_path("gc-same-q")),
        render_double_pendulum_figure(double, figure_path("gc-double-pendulum-config")),
    ]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.add_argument("--theta1-deg", type=float, default=90.0, help="二重振り子の初期角 θ1 [deg]")
    p.add_argument("--theta2-deg", type=float, default=0.0, help="二重振り子の初期角 θ2 [deg]")
    p.add_argument("--duration", type=float, default=10.0, help="二重振り子の長さ [s]")
    p.add_argument("--fps", type=int, default=30, help="動画のフレームレート")
    p.add_argument("--no-video", action="store_true", help="動画を作らず、図だけ作る")
    args = p.parse_args(argv)

    runs = same_position_different_velocity(np.deg2rad(60.0), np.array([0.0, 2.0, -3.0, 4.0]))
    print(f"[{NAME}] 単振り子 θ(0)=60° で θ̇(0) だけ変える")
    for w, measured, theory in zip(
        runs.omega0s, runs.turning_angle_measured, runs.turning_angle_energy, strict=True
    ):
        print(
            f"  θ̇(0)={w:+.0f} rad/s: 折り返し 実測 {np.rad2deg(measured):.3f}°"
            f" / エネルギー保存 {np.rad2deg(theory):.3f}°"
        )

    params = DoublePendulumParams(
        theta1=np.deg2rad(args.theta1_deg), theta2=np.deg2rad(args.theta2_deg)
    )
    run = simulate_double_pendulum(params, args.duration)
    drift = np.max(np.abs(run.energy - run.energy[0])) / abs(run.energy[0])
    print(f"二重振り子: エネルギーの最大相対誤差 {drift:.1e}（RK4, dt=1 ms, {args.duration:g} s）")
    for mode in (0, 1):
        measured, theory = measure_normal_mode_period(DoublePendulumParams(), mode)
        print(f"  基準振動 {mode + 1}: 周期 実測 {measured:.5f} s / 理論 {theory:.5f} s")

    print("図を作成中...")
    for path in figures():  # 引数を変えても、教科書の図は既定の条件のまま
        print(f"  -> {path}")
    if not args.no_video:
        print("動画を作成中...")
        video = render_double_pendulum_video(
            run, output_path(TOPIC, "double_pendulum", ".mp4"), args.fps
        )
        print(f"  -> {video}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
