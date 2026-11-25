"""`physlab run areal-velocity`: 中心力の軌道は、同じ時間に同じ面積を掃く。

設計書: docs/design/006-conservation-laws.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from physlab.application.mechanics.areal_velocity import (
    central_orbit,
    equal_time_sectors,
    time_of_one_revolution,
)
from physlab.domain.mechanics.central import CentralParams
from physlab.infrastructure.rendering.mechanics.areal_velocity import (
    render_areal_figure,
    render_areal_video,
)
from physlab.infrastructure.storage import figure_path, output_path

NAME = "areal-velocity"
TITLE = "角運動量と面積速度：中心力なら、同じ時間に掃く面積は同じ"
TOPIC = "mechanics"


def _sectors(exponent: float):
    orbit = central_orbit(CentralParams(exponent=exponent))
    # 扇形が重ならないように、中心を 1 周するまでの時間を 12 等分する
    return equal_time_sectors(orbit, 12, span=time_of_one_revolution(orbit))


def figures() -> list[Path]:
    """教科書の図 1 枚（``outputs/figures/am-areal.pdf``）。"""
    return [render_areal_figure(_sectors(CentralParams().exponent), figure_path("am-areal"))]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.add_argument("--exponent", type=float, default=1.4, help="U = −1/rⁿ の n（0 < n < 2）")
    p.add_argument("--fps", type=int, default=30, help="動画のフレームレート")
    p.add_argument("--no-video", action="store_true", help="動画を作らず、図だけ作る")
    args = p.parse_args(argv)

    sectors = _sectors(args.exponent)
    orbit = sectors.orbit
    radius = np.linalg.norm(orbit.q, axis=1)
    print(f"[{NAME}] U = −1/r^{args.exponent:g}、r の範囲 {radius.min():.3f} – {radius.max():.3f}")
    low, high = sectors.areas.min(), sectors.areas.max()
    print(f"  1 区間 {sectors.duration:.4f} s の面積: {low:.10f} – {high:.10f}")
    print(f"  弧の長さ: {sectors.arc_lengths.min():.3f} – {sectors.arc_lengths.max():.3f}")
    print(f"  角運動量のずれの最大 = {np.max(np.abs(orbit.angular - orbit.angular[0])):.1e}")

    print("図を作成中...")
    for path in figures():  # 引数を変えても、教科書の図は既定の条件のまま
        print(f"  -> {path}")
    if not args.no_video:
        print("動画を作成中...")
        video = render_areal_video(
            sectors, output_path(TOPIC, "areal_velocity", ".mp4"), fps=args.fps
        )
        print(f"  -> {video}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
