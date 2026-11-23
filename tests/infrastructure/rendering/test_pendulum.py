"""描画の煙テスト: ファイルができることだけ確認する（見た目は動画のフレームを目視で確認する）。"""

import numpy as np

from physlab.application.mechanics.pendulum import compare_energy_drift, compare_with_small_angle
from physlab.domain.mechanics.pendulum import PendulumParams
from physlab.infrastructure.rendering.mechanics.pendulum import (
    render_comparison_video,
    render_energy_drift_figure,
)


def test_render_comparison_video_writes_file(tmp_path) -> None:
    params = PendulumParams(theta0=np.deg2rad(120.0))
    result = compare_with_small_angle(params, duration=0.5)
    saved = render_comparison_video(result, tmp_path / "v.mp4", fps=10)
    assert saved.suffix in {".mp4", ".gif"}  # ffmpeg が無い環境では GIF になる
    assert saved.stat().st_size > 0


def test_render_energy_drift_figure_writes_file(tmp_path) -> None:
    params = PendulumParams(theta0=np.deg2rad(120.0))
    result = compare_energy_drift(params, duration=2.0, dt=0.01)
    saved = render_energy_drift_figure(result, tmp_path / "e.png")
    assert saved.stat().st_size > 0
