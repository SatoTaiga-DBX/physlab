"""`physlab run pendulum`: 単振り子。引数を解釈して application と infrastructure をつなぐ。

設計書: docs/design/002-pendulum.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from physlab.application.mechanics.pendulum import (
    compare_energy_drift,
    compare_with_small_angle,
    period_ratio_curve,
    phase_portrait,
)
from physlab.domain.mechanics.pendulum import PendulumParams
from physlab.infrastructure.rendering.mechanics.pendulum import (
    render_comparison_video,
    render_energy_drift_figure,
    render_period_ratio_figure,
    render_phase_portrait_figure,
)
from physlab.infrastructure.storage import figure_path, output_path

NAME = "pendulum"
TITLE = "単振り子：大振幅での小角近似の破綻と、積分法のエネルギー保存"
TOPIC = "mechanics"


def figures() -> list[Path]:
    """教科書の図 3 枚（``outputs/figures/pendulum-*.pdf``）。"""
    energy = compare_energy_drift(PendulumParams(), 100.0, 0.01)
    return [
        render_period_ratio_figure(
            # 注記する 30°・120° を格子点に含める（近い点の値でなく、ちょうどの値を示すため）
            period_ratio_curve(np.union1d(np.linspace(1.0, 179.0, 300), [30.0, 120.0])),
            figure_path("pendulum-period-ratio"),
        ),
        render_phase_portrait_figure(phase_portrait(), figure_path("pendulum-phase-portrait")),
        render_energy_drift_figure(energy, figure_path("pendulum-energy-drift")),
    ]


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.add_argument("--theta0-deg", type=float, default=120.0, help="初期角 [deg]（既定 120）")
    p.add_argument("--length", type=float, default=1.0, help="振り子の長さ [m]")
    p.add_argument("--duration", type=float, default=10.0, help="動画の長さ [s]")
    p.add_argument("--fps", type=int, default=30, help="動画のフレームレート")
    p.add_argument("--energy-duration", type=float, default=100.0, help="エネルギー比較の長さ [s]")
    p.add_argument("--energy-dt", type=float, default=0.01, help="エネルギー比較の刻み幅 [s]")
    args = p.parse_args(argv)

    params = PendulumParams(length=args.length, theta0=np.deg2rad(args.theta0_deg))
    print(f"[{NAME}] θ0={args.theta0_deg}°, L={args.length} m")

    print("動画を作成中...")
    comparison = compare_with_small_angle(params, args.duration)
    print(
        f"  周期 厳密 {comparison.period_exact:.4f} s / 小角近似 {comparison.period_linear:.4f} s"
    )
    video = render_comparison_video(comparison, output_path(TOPIC, NAME, ".mp4"), fps=args.fps)
    print(f"  -> {video}")

    print("エネルギー保存の比較図を作成中...")
    energy = compare_energy_drift(params, args.energy_duration, args.energy_dt)
    figure = render_energy_drift_figure(energy, output_path(TOPIC, f"{NAME}_energy_drift", ".png"))
    print(f"  -> {figure}")
    for method, err in energy.relative_error.items():
        print(f"  {method:>6}: 最終 {err[-1]:+.2e} / 最大絶対誤差 {np.abs(err).max():.2e}")


if __name__ == "__main__":
    main()
