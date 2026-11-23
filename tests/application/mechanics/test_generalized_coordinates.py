"""一般座標: 位置だけでは運動が決まらないこと、二重振り子のエネルギー保存と基準振動。"""

import numpy as np
import pytest

from physlab.application.mechanics.generalized_coordinates import (
    measure_normal_mode_period,
    same_position_different_velocity,
    simulate_double_pendulum,
)
from physlab.domain.mechanics.double_pendulum import DoublePendulumParams


def test_same_position_different_velocity_gives_different_motions():
    theta0 = np.deg2rad(60.0)
    omega0s = np.array([0.0, 2.0, -3.0, 4.0])
    runs = same_position_different_velocity(theta0, omega0s, duration=6.0)
    assert runs.theta[:, 0] == pytest.approx(theta0)  # 出発点は同じ
    assert runs.omega[:, 0] == pytest.approx(omega0s)
    # 1 秒後にはみな違う位置にいる
    at_1s = runs.theta[:, np.searchsorted(runs.t, 1.0)]
    assert np.min(np.abs(np.subtract.outer(at_1s, at_1s) + np.eye(4))) > 0.1


def test_turning_angle_follows_energy_conservation():
    """cos θ_turn = cos θ0 − L ω0²/(2g)。符号の違う ω0 では折り返しの角度が同じ。"""
    runs = same_position_different_velocity(np.deg2rad(60.0), np.array([0.0, 2.0, -2.0, 4.0]))
    assert runs.turning_angle_measured == pytest.approx(runs.turning_angle_energy, rel=1e-4)
    assert runs.turning_angle_measured[1] == pytest.approx(runs.turning_angle_measured[2], rel=1e-6)
    assert np.rad2deg(runs.turning_angle_energy[0]) == pytest.approx(60.0)


def test_fast_start_goes_over_the_top():
    """E/m > 2gL なら回転運動になり、折り返しの角度は存在しない（NaN）。"""
    runs = same_position_different_velocity(np.deg2rad(60.0), np.array([6.0]), duration=3.0)
    assert np.isnan(runs.turning_angle_energy[0])
    assert runs.turning_angle_measured[0] > np.pi  # 実際に真上を越えている


def test_double_pendulum_conserves_energy():
    run = simulate_double_pendulum(DoublePendulumParams(), duration=10.0, dt=1e-3)
    relative = np.abs(run.energy - run.energy[0]) / abs(run.energy[0])
    assert relative.max() < 1e-7
    # 動いている（エネルギーが保存されるだけの静止ではない）
    assert np.ptp(run.theta2) > 1.0


@pytest.mark.parametrize("mode", [0, 1])
def test_small_oscillation_period_matches_normal_mode_theory(mode: int):
    measured, theory = measure_normal_mode_period(DoublePendulumParams(), mode)
    assert measured == pytest.approx(theory, rel=1e-3)
