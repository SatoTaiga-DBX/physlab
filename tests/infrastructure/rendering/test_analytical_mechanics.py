"""解析力学の動画の煙テスト（短く小さい設定で、動画ファイルが書き出されること）。"""

from pathlib import Path

import numpy as np

from physlab.application.mechanics.galilean import simulate_two_body_frames
from physlab.application.mechanics.generalized_coordinates import simulate_double_pendulum
from physlab.application.mechanics.least_action import relax_path
from physlab.domain.mechanics.double_pendulum import DoublePendulumParams
from physlab.domain.mechanics.lagrangian import uniform_gravity
from physlab.domain.mechanics.two_body import TwoBodyParams
from physlab.infrastructure.rendering.mechanics.galilean import render_two_body_video
from physlab.infrastructure.rendering.mechanics.generalized_coordinates import (
    render_double_pendulum_video,
)
from physlab.infrastructure.rendering.mechanics.least_action import render_relaxation_video

VIDEO_SUFFIXES = {".mp4", ".gif"}


def test_relaxation_video_is_written(tmp_path: Path) -> None:
    relax = relax_path(
        uniform_gravity(), 0.0, 0.0, 2.0, 0.0, n_iter=40, n_frames=6, frame_horizon=40
    )
    path = render_relaxation_video(relax, tmp_path / "relax.mp4", fps=5, repeat=1, hold_seconds=0.2)
    assert path.suffix in VIDEO_SUFFIXES and path.stat().st_size > 0


def test_double_pendulum_video_is_written(tmp_path: Path) -> None:
    run = simulate_double_pendulum(DoublePendulumParams(), duration=0.6, dt=1e-3)
    path = render_double_pendulum_video(run, tmp_path / "dp.mp4", fps=5)
    assert path.suffix in VIDEO_SUFFIXES and path.stat().st_size > 0


def test_two_body_video_is_written(tmp_path: Path) -> None:
    frames = simulate_two_body_frames(TwoBodyParams(), np.array([0.3, 0.1]), 0.6, dt=1e-3)
    path = render_two_body_video(frames, tmp_path / "tb.mp4", fps=5)
    assert path.suffix in VIDEO_SUFFIXES and path.stat().st_size > 0
