"""単振り子の物理（domain）: 時間積分を使わずに、法則そのもの（式・極限・特殊値）を確かめる。"""

import numpy as np
import pytest

from physlab.domain.mechanics.pendulum import (
    PendulumParams,
    deriv,
    energy,
    period_exact,
    period_small_angle,
    small_angle_solution,
)


def test_small_angle_period_formula() -> None:
    params = PendulumParams(g=9.81, length=1.0)
    assert period_small_angle(params) == pytest.approx(2 * np.pi * np.sqrt(1.0 / 9.81))
    assert period_small_angle(params) == pytest.approx(2.0061, abs=1e-4)  # 教科書の値


def test_exact_period_tends_to_small_angle_period() -> None:
    params = PendulumParams(theta0=np.deg2rad(1.0))
    assert period_exact(params) == pytest.approx(period_small_angle(params), rel=1e-4)


def test_exact_period_known_value_at_120_degrees() -> None:
    params = PendulumParams(theta0=np.deg2rad(120.0))
    assert period_exact(params) / period_small_angle(params) == pytest.approx(1.3729, abs=1e-4)


def test_exact_period_grows_with_amplitude() -> None:
    periods = [period_exact(PendulumParams(theta0=np.deg2rad(a))) for a in (10, 60, 120, 170)]
    assert periods == sorted(periods)


def test_exact_period_series_expansion() -> None:
    """T ≈ T0 (1 + θ0²/16) が小さい振幅でよく成り立つ。"""
    params = PendulumParams(theta0=np.deg2rad(20.0))
    approx = period_small_angle(params) * (1 + params.theta0**2 / 16)
    assert period_exact(params) == pytest.approx(approx, rel=1e-4)


def test_energy_at_release_and_at_bottom() -> None:
    params = PendulumParams(g=9.81, length=2.0, theta0=np.deg2rad(60.0))
    release = energy(np.array(params.theta0), np.array(0.0), params)
    assert release == pytest.approx(9.81 * 2.0 * (1 - np.cos(params.theta0)))
    # 最下点では全て運動エネルギー: ½ (L ω)²
    omega_bottom = np.sqrt(2 * release) / params.length
    assert energy(np.array(0.0), np.array(omega_bottom), params) == pytest.approx(release)


def test_equation_of_motion_special_points() -> None:
    params = PendulumParams(g=9.81, length=1.0)
    f = deriv(params)
    assert f(0.0, np.array([0.0, 2.0])) == pytest.approx([2.0, 0.0])  # 最下点: 加速度ゼロ
    assert f(0.0, np.array([np.pi / 2, 0.0])) == pytest.approx([0.0, -9.81])  # 水平: -g/L
    assert f(0.0, np.array([np.pi, 0.0])) == pytest.approx([0.0, 0.0], abs=1e-12)  # 真上: つり合い


def test_small_angle_solution_is_periodic() -> None:
    params = PendulumParams(theta0=np.deg2rad(5.0))
    t = np.array([0.0, period_small_angle(params)])
    theta = small_angle_solution(t, params)
    assert theta[0] == pytest.approx(params.theta0)
    assert theta[1] == pytest.approx(params.theta0)
