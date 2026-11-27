"""中心力ポテンシャル: 微分の式と、式の切り替わる点の連続性。"""

import numpy as np
import pytest

from physlab.domain.mechanics.potentials import (
    Gaussian,
    PowerLaw,
    SoftDisk,
    UniformBall,
    Yukawa,
)

POTENTIALS = [
    PowerLaw(alpha=0.5, n=1.0),
    PowerLaw(alpha=-0.5, n=1.0),
    PowerLaw(alpha=0.1, n=12.0),
    Gaussian(height=1.0, width=0.8),
    UniformBall(alpha=0.5, radius=1.4),
    Yukawa(alpha=1.0, length=0.6),
    SoftDisk(stiffness=4.0, radius=1.0),
]


@pytest.mark.parametrize("potential", POTENTIALS, ids=lambda p: type(p).__name__)
def test_derivative_matches_a_numerical_difference(potential) -> None:
    r = np.array([0.7, 0.95, 1.2, 1.7, 3.0])
    h = 1e-6
    numerical = (potential.value(r + h) - potential.value(r - h)) / (2 * h)
    assert potential.derivative(r) == pytest.approx(numerical, rel=1e-6, abs=1e-8)


def test_uniform_ball_joins_the_point_charge_smoothly() -> None:
    ball = UniformBall(alpha=0.5, radius=2.0)
    edge = np.array([2.0 - 1e-9, 2.0 + 1e-9])
    assert ball.value(edge) == pytest.approx(0.25, abs=1e-8)
    assert ball.derivative(edge) == pytest.approx(-0.125, abs=1e-8)
    assert float(ball.value(np.array(0.0))) == pytest.approx(1.5 * 0.5 / 2.0)  # 中心は表面の 1.5 倍
    assert ball.breakpoint_radius == 2.0
