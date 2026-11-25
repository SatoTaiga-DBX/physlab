"""`physlab run symmetry`: 対称性と保存量（外場で対称性を壊す、積分法、回転する棒）。

設計書: docs/design/006-conservation-laws.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from physlab.application.mechanics.symmetry import (
    compare_integrators,
    rotating_rod_run,
    symmetry_experiment,
)
from physlab.infrastructure.rendering.mechanics.symmetry import (
    render_cluster_video,
    render_integrator_figure,
    render_rod_figure,
    render_symmetry_figure,
)
from physlab.infrastructure.storage import figure_path, output_path

NAME = "symmetry"
TITLE = "対称性と保存量：外場で対称性を壊すと、どの量が保存しなくなるか"
TOPIC = "mechanics"


def figures() -> list[Path]:
    """教科書の図 3 枚（``outputs/figures/cons-*.pdf``）。常に教科書の条件で作る。"""
    return [
        render_symmetry_figure(symmetry_experiment(10.0), figure_path("cons-symmetry")),
        render_integrator_figure(compare_integrators(), figure_path("cons-integrators")),
        render_rod_figure(rotating_rod_run(), figure_path("cons-rotating-rod")),
    ]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.add_argument("--duration", type=float, default=10.0, help="実験の時間の長さ")
    p.add_argument("--fps", type=int, default=30, help="動画のフレームレート")
    p.add_argument("--no-video", action="store_true", help="動画を作らず、図だけ作る")
    args = p.parse_args(argv)

    runs = symmetry_experiment(args.duration)
    print(f"[{NAME}] 3 質点をばねでつないだ系の、保存量の変化の最大値（t = 0 からのずれ）")
    print("  外場            |ΔE|      |ΔPx|     |ΔPy|     |ΔL|")
    for field, run in runs.items():
        values = (run.energy, run.momentum[:, 0], run.momentum[:, 1], run.angular)
        changes = "  ".join(f"{np.max(np.abs(v - v[0])):.1e}" for v in values)
        print(f"  {field:<13} {changes}")

    c = compare_integrators()
    print(f"積分法の比較（外場なし、刻み {c.dt}、時間 {c.t[-1]:g}）の最大のずれ")
    for method in ("verlet", "rk4"):
        print(
            f"  {method:<6} E {c.energy_error[method].max():.1e}"
            f"  P {c.momentum_error[method].max():.1e}  L {c.angular_error[method].max():.1e}"
        )

    rod = rotating_rod_run()
    e = rod.energy_function
    print(f"回転する棒: E = pṙ − L は {e[0]:.6f} → {e[-1]:.6f}（保存）")
    print(
        f"  T は {rod.kinetic[0]:.4f} → {rod.kinetic[-1]:.4f}、"
        f"モーターの仕事 {rod.motor_work[-1]:.4f}"
    )

    print("図を作成中...")
    for path in figures():  # 引数を変えても、教科書の図は既定の条件のまま
        print(f"  -> {path}")
    if not args.no_video:
        print("動画を作成中...")
        video = render_cluster_video(
            runs["none"], output_path(TOPIC, "symmetry_cluster", ".mp4"), fps=args.fps
        )
        print(f"  -> {video}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
