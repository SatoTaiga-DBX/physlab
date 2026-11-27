"""`physlab run decay`: 粒子の崩壊。破片の向きとエネルギーの分布を保存則から。

設計書: docs/design/007-collisions.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from physlab.application.mechanics.decay import (
    angle_histogram,
    decay_burst,
    energy_error,
    energy_histogram,
    momentum_error,
    simulate_decays,
)
from physlab.domain.mechanics.decay import DecayParams, kinetic_energy_range, max_lab_angle
from physlab.infrastructure.rendering.mechanics.decay import (
    render_burst_video,
    render_distribution_figure,
    render_velocity_figure,
)
from physlab.infrastructure.storage import figure_path, output_path

NAME = "decay"
TITLE = "粒子の崩壊：破片の向きとエネルギーの分布は、保存則で決まる"
TOPIC = "mechanics"
RATIOS = (0.5, 2.0)  # 親の速さ / 静止系での破片の速さ


def _samples(n: int, seed: int):
    params = DecayParams()
    v10 = params.rest_speeds[0]
    return [simulate_decays(params, k * v10, n=n, seed=seed) for k in RATIOS]


def figures() -> list[Path]:
    """教科書の図 2 枚。"""
    cases = [
        (k, angle_histogram(s), energy_histogram(s))
        for k, s in zip(RATIOS, _samples(200_000, 1), strict=True)
    ]
    return [
        render_velocity_figure(RATIOS, figure_path("dec-velocity")),
        render_distribution_figure(cases, figure_path("dec-distributions")),
    ]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.add_argument("--samples", type=int, default=200_000, help="崩壊させる回数")
    p.add_argument("--seed", type=int, default=1, help="乱数の種")
    p.add_argument("--fps", type=int, default=30, help="動画のフレームレート")
    p.add_argument("--no-video", action="store_true", help="動画を作らず、図だけ作る")
    args = p.parse_args(argv)

    params = DecayParams()
    v10, v20 = params.rest_speeds
    print(f"[{NAME}] m1={params.m1:g} kg, m2={params.m2:g} kg, 解放エネルギー {params.energy:g} J")
    print(f"  静止系: 運動量 {params.momentum:g} kg m/s、速さ {v10:g} と {v20:g} m/s")
    for k, sample in zip(RATIOS, _samples(args.samples, args.seed), strict=True):
        angle = angle_histogram(sample)
        energy = energy_histogram(sample)
        low, high = kinetic_energy_range(params.m1, sample.v_parent, v10)
        print(f"  V = {sample.v_parent:g} m/s（V/v0 = {k:g}）:")
        print(
            f"    角の上限 {np.degrees(max_lab_angle(k)):.1f}°、"
            f"実測の最大 {np.degrees(sample.theta.max()):.2f}°、"
            f"前方（90°未満）の割合 {np.mean(sample.theta < np.pi / 2):.4f}"
        )
        print(f"    運動エネルギー {low:g} – {high:g} J、平均 {sample.kinetic.mean():.3f} J")
        z_angle = np.max(np.abs(angle.measured - angle.expected) / angle.error)
        z_energy = np.max(np.abs(energy.measured - energy.expected) / energy.error)
        print(f"    理論との最大のずれ（統計誤差の倍数）: 角 {z_angle:.1f}、E {z_energy:.1f}")
        print(
            f"    保存則のずれ: 運動量 {momentum_error(sample):.1e}、"
            f"エネルギー {energy_error(sample):.1e}"
        )

    print("図を作成中...")
    for path in figures():
        print(f"  -> {path}")
    if not args.no_video:
        print("動画を作成中...")
        unit = DecayParams(energy=2.0 / 3.0)  # 静止系の破片の速さが 1 になる
        t = np.arange(-1.0, 3.0 + 0.005, 0.01)
        bursts = [decay_burst(unit, k * unit.rest_speeds[0], t) for k in RATIOS]
        video = render_burst_video(bursts, output_path(TOPIC, "decay_cone", ".mp4"), fps=args.fps)
        print(f"  -> {video}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
