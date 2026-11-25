"""`physlab run center-of-mass`: 回しながら投げた 3 質点の重心は放物線を描く。

設計書: docs/design/006-conservation-laws.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from physlab.application.mechanics.center_of_mass import throw_cluster
from physlab.infrastructure.rendering.mechanics.center_of_mass import (
    render_decomposition_figure,
    render_throw_figure,
    render_throw_video,
)
from physlab.infrastructure.storage import figure_path, output_path

NAME = "center-of-mass"
TITLE = "慣性中心：回りながら振動して飛んでも、重心は放物線を描く"
TOPIC = "mechanics"


def figures() -> list[Path]:
    """教科書の図 2 枚（``outputs/figures/com-*.pdf``）。"""
    result = throw_cluster()
    return [
        render_throw_figure(result, figure_path("com-throw")),
        render_decomposition_figure(result, figure_path("com-decomposition")),
    ]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.add_argument("--spin", type=float, default=6.0, help="投げるときの回転の角速度 [rad/s]")
    p.add_argument("--fps", type=int, default=30, help="動画のフレームレート")
    p.add_argument("--no-video", action="store_true", help="動画を作らず、図だけ作る")
    args = p.parse_args(argv)

    r = throw_cluster(spin=args.spin)
    print(f"[{NAME}] 初速 (3, 8) m/s、回転 {args.spin:g} rad/s で投げる（飛行 {r.t[-1]:.3f} s）")
    print(f"  重心と放物線の式の差の最大 = {np.max(np.abs(r.center - r.center_exact)):.1e} m")
    internal = r.internal_energy
    drift = np.max(np.abs(internal / internal[0] - 1))
    print(f"  内部エネルギー {internal[0]:.4f} J の相対変化の最大 = {drift:.1e}")
    spin = r.angular_relative
    print(
        f"  重心まわりの角運動量 {spin[0]:.4f} のずれの最大 = {np.max(np.abs(spin - spin[0])):.1e}"
    )
    print(
        f"  原点まわりの角運動量 {r.angular_total[0]:.3f} → {r.angular_total[-1]:.3f}（保存しない）"
    )

    print("図を作成中...")
    for path in figures():
        print(f"  -> {path}")
    if not args.no_video:
        print("動画を作成中...")
        video = render_throw_video(r, output_path(TOPIC, "com_throw", ".mp4"), fps=args.fps)
        print(f"  -> {video}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
