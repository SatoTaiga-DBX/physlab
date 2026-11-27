"""衝突・散乱の動画の煙テスト（短く小さい設定で、動画ファイルが書き出されること）。"""

from pathlib import Path

import numpy as np

from physlab.application.mechanics.decay import decay_burst
from physlab.application.mechanics.deflection import deflection, trace_orbits
from physlab.application.mechanics.elastic_collision import (
    DEMO_POTENTIALS,
    simulate_collisions,
)
from physlab.application.mechanics.scattering import deflection_curve, hill
from physlab.domain.mechanics.decay import DecayParams
from physlab.domain.mechanics.potentials import coulomb
from physlab.infrastructure.rendering.mechanics.decay import render_burst_video
from physlab.infrastructure.rendering.mechanics.elastic_collision import render_collision_video
from physlab.infrastructure.rendering.mechanics.rutherford import render_rutherford_video
from physlab.infrastructure.rendering.mechanics.scattering import render_beam_video

VIDEO_SUFFIXES = {".mp4", ".gif"}


def _written(path: Path) -> bool:
    return path.suffix in VIDEO_SUFFIXES and path.stat().st_size > 0


def test_decay_video_is_written(tmp_path: Path) -> None:
    t = np.arange(-0.2, 0.41, 0.01)
    bursts = [decay_burst(DecayParams(), v, t, n=12) for v in (3.0, 12.0)]
    assert _written(render_burst_video(bursts, tmp_path / "decay.mp4", fps=5))


def test_collision_video_is_written(tmp_path: Path) -> None:
    batch = simulate_collisions(
        1.0, 2.0, DEMO_POTENTIALS["gaussian"], [0.5], start=1.0, duration=1.0, dt=5e-3
    )
    assert _written(render_collision_video(batch, tmp_path / "col.mp4", fps=5))


def test_beam_video_is_written(tmp_path: Path) -> None:
    rho = np.array([-1.0, -0.3, 0.3, 1.0])
    orbits = trace_orbits(hill(), 1.0, 2.0, rho, start_radius=2.0, duration=1.0, dt=5e-3)
    chi = np.abs(deflection(np.abs(rho), hill(), 2.0))
    curve = deflection_curve(hill(), 2.0, 1.5, n=50)
    assert _written(render_beam_video(orbits, chi, curve, tmp_path / "beam.mp4", fps=5))


def test_rutherford_video_is_written(tmp_path: Path) -> None:
    orbits = trace_orbits(
        coulomb(0.5), 1.0, 1.0, [-1.0, 1.0], start_radius=3.0, duration=1.0, dt=5e-3
    )
    assert _written(render_rutherford_video(orbits, 1.0, tmp_path / "ruth.mp4", fps=5))
