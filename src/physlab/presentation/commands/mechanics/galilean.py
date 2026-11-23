"""`physlab run galilean`: ガリレイの相対性（電車の投げ上げ、作用のずれ、二体問題）。

設計書: docs/design/005-analytical-mechanics-ch1.md
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from physlab.application.mechanics.galilean import (
    action_shift_under_boost,
    simulate_two_body_frames,
    throw_in_train,
)
from physlab.domain.mechanics.two_body import TwoBodyParams
from physlab.infrastructure.rendering.mechanics.galilean import (
    render_action_shift_figure,
    render_train_figure,
    render_two_body_figure,
    render_two_body_video,
)
from physlab.infrastructure.storage import figure_path, output_path

NAME = "galilean"
TITLE = "ガリレイの相対性：運動方程式は形を変えず、作用は定数だけずれる（電車・二体問題）"
TOPIC = "mechanics"

DEFAULT_DRIFT = (0.35, 0.2)  # 重心系 K' の実験室系 K に対する速度


def _two_body(drift: tuple[float, float], duration: float | None):
    params = TwoBodyParams()
    return simulate_two_body_frames(
        params, np.array(drift), 1.5 * params.period if duration is None else duration
    )


def figures() -> list[Path]:
    """教科書の図 3 枚（``outputs/figures/gal-*.pdf``）。"""
    train = throw_in_train(train_speed=10.0, throw_speed=7.0)
    shift = action_shift_under_boost(2.0, 3.0, 0.0, 0.0, 2.0, 1.0, np.linspace(-1.5, 1.5, 13))
    return [
        render_train_figure(train, figure_path("gal-train")),
        render_action_shift_figure(shift, figure_path("gal-action-shift")),
        render_two_body_figure(_two_body(DEFAULT_DRIFT, None), figure_path("gal-two-body")),
    ]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog=f"physlab run {NAME}", description=TITLE)
    p.add_argument("--drift-x", type=float, default=DEFAULT_DRIFT[0], help="K' の速度の x 成分")
    p.add_argument("--drift-y", type=float, default=DEFAULT_DRIFT[1], help="K' の速度の y 成分")
    p.add_argument(
        "--duration", type=float, default=None, help="二体問題の長さ（既定は周期の 1.5 倍）"
    )
    p.add_argument("--fps", type=int, default=30, help="動画のフレームレート")
    p.add_argument("--no-video", action="store_true", help="動画を作らず、図だけ作る")
    args = p.parse_args(argv)

    train = throw_in_train(10.0, 7.0)
    print(f"[{NAME}] 電車（V=10 m/s）の中で 7 m/s で真上に投げる")
    print(f"  加速度 K': {train.accel_train[0]}, K: {train.accel_ground[0]}（どちらも (0, −g)）")
    dke = train.kinetic_ground[0] - train.kinetic_train[0]
    print(f"  運動エネルギーの差 K−K' = {dke:.3f}（= ½V² = 50）")

    shift = action_shift_under_boost(2.0, 3.0, 0.0, 0.0, 2.0, 1.0, np.linspace(-1.5, 1.5, 13))
    lo, hi = shift.shift.min(), shift.shift.max()
    print(f"  作用のずれ S−S' = {lo:.6f} … {hi:.6f}（理論 {shift.shift_theory:.6f}）")

    drift = (args.drift_x, args.drift_y)
    frames = _two_body(drift, args.duration)
    params = frames.params
    print(f"二体問題 m1={params.m1:g}, m2={params.m2:g}, 周期 T={params.period:.3f}, V={drift}")
    print(
        f"  K で独立に解いた軌道 − (K' の軌道 + V t) の最大 = {frames.frame_difference.max():.1e}"
    )
    e0 = frames.energy_cm[0]
    print(
        f"  エネルギーの最大相対誤差（K'）= {np.max(np.abs(frames.energy_cm - e0)) / abs(e0):.1e}"
    )
    print(
        f"  K の全運動量 = {frames.momentum_lab[-1]}（M V = {params.total_mass * np.array(drift)}）"
    )

    print("図を作成中...")
    for path in figures():
        print(f"  -> {path}")
    if not args.no_video:
        print("動画を作成中...")
        video = render_two_body_video(
            frames, output_path(TOPIC, "galilean_two_body", ".mp4"), args.fps
        )
        print(f"  -> {video}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
