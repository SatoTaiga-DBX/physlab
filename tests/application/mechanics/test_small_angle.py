"""力積の近似: 公式との一致と、誤差の減り方。"""

import numpy as np
import pytest

from physlab.application.mechanics.small_angle import (
    FAMILIES,
    error_sweep,
    impulse_angle,
    impulse_demo,
    small_angle_tail,
)
from physlab.domain.mechanics.potentials import Gaussian, PowerLaw, UniformBall
from physlab.domain.mechanics.scattering import ball_small_angle, power_law_small_angle

RHO = np.array([0.4, 1.0, 3.0])


@pytest.mark.parametrize("n", [1.0, 2.0, 4.0])
def test_impulse_integral_matches_the_gamma_function_formula(n: float) -> None:
    theta = impulse_angle(RHO, PowerLaw(0.5, n))
    assert theta == pytest.approx(power_law_small_angle(RHO, 0.5, n, 1.0, 1.0), rel=1e-10)


def test_impulse_integral_for_the_gaussian_hill() -> None:
    height, width, energy = 0.02, 1.0, 0.5
    theta = impulse_angle(RHO, Gaussian(height, width))
    expected = np.sqrt(np.pi) * height / energy * RHO / width * np.exp(-((RHO / width) ** 2))
    assert theta == pytest.approx(expected, rel=1e-9)


def test_impulse_integral_for_the_ball_is_the_projected_charge() -> None:
    rho = np.array([0.3, 1.0, 1.9, 2.5])
    theta = impulse_angle(rho, UniformBall(0.5, 2.0))
    assert theta == pytest.approx(ball_small_angle(rho, 1.0, 2.0), rel=1e-6)


@pytest.mark.parametrize("family", list(FAMILIES))
def test_error_decreases_with_the_expected_power(family: str) -> None:
    sweep = error_sweep(family)
    assert sweep.relative_error[0] < 0.02
    # クーロン力だけは引力と斥力が鏡像になるので、1 次の誤差が消える
    expected = 2.0 if family == "coulomb" else 1.0
    assert sweep.slope == pytest.approx(expected, abs=0.05)


@pytest.mark.parametrize("n", [1.0, 2.0, 4.0])
def test_small_angle_cross_section_exponent(n: float) -> None:
    tail = small_angle_tail(n)
    assert tail.slope == pytest.approx(-2.0 - 2.0 / n, abs=0.02)


def test_demo_orbit_is_nearly_straight() -> None:
    demo = impulse_demo()
    assert demo.approx_angle == pytest.approx(demo.exact_angle, rel=0.03)
    # 横向きの運動量の変化は、力を時間積分したもの
    impulse = np.trapezoid(demo.force_y, demo.t)
    assert demo.momentum_y[-1] - demo.momentum_y[0] == pytest.approx(impulse, abs=1e-7)
    straight = np.trapezoid(demo.straight_force_y, demo.straight_x) / np.sqrt(20.0)
    assert demo.approx_momentum_y == pytest.approx(straight, rel=1e-6)
