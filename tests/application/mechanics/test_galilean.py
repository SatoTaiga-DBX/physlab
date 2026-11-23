"""ガリレイの相対性: 運動方程式の不変性、作用のずれ、二体問題の 2 つの慣性系。"""

import numpy as np
import pytest

from physlab.application.mechanics.galilean import (
    action_shift_under_boost,
    simulate_two_body_frames,
    throw_in_train,
)
from physlab.domain.mechanics.two_body import TwoBodyParams


def test_train_throw_is_the_same_motion_in_both_frames():
    r = throw_in_train(train_speed=10.0, throw_speed=7.0)
    # K の軌跡 = K' の軌跡 + V t
    expected = r.train + np.stack([10.0 * r.t, np.zeros_like(r.t)], axis=1)
    assert r.ground == pytest.approx(expected, abs=1e-9)
    # どちらの系でも加速度は (0, −g)（運動方程式は形を変えない）
    for a in (r.accel_train, r.accel_ground):
        assert a[:, 0] == pytest.approx(0.0, abs=1e-6)
        assert a[:, 1] == pytest.approx(-9.81, abs=1e-6)
    # K' では真上に上がって降りる（解析解 y = v0 t − g t²/2 と一致）。K では同じ時間に V t だけ進む
    t_end = r.t[-1]  # 時間格子の最終点（滞空時間 2 v0/g にほぼ等しい）
    assert t_end == pytest.approx(2 * 7.0 / 9.81, abs=1e-3)
    assert r.train[:, 0] == pytest.approx(0.0, abs=1e-12)
    assert r.train[-1, 1] == pytest.approx(7.0 * t_end - 0.5 * 9.81 * t_end**2, abs=1e-9)
    assert abs(r.train[-1, 1]) < 1e-2  # 投げ上げた高さに戻っている
    assert r.ground[-1, 0] == pytest.approx(10.0 * t_end, rel=1e-9)


def test_kinetic_energy_is_not_invariant_but_the_difference_is_constant():
    r = throw_in_train(train_speed=10.0, throw_speed=7.0)
    assert not np.allclose(r.kinetic_ground, r.kinetic_train)
    assert r.kinetic_ground - r.kinetic_train == pytest.approx(0.5 * 10.0**2, rel=1e-9)


def test_action_shift_is_independent_of_the_path():
    eps = np.linspace(-1.5, 1.5, 13)
    s = action_shift_under_boost(2.0, 3.0, 0.0, 0.0, 2.0, 1.0, eps)
    assert s.shift_theory == pytest.approx(-12.0)  # m V Δq' + ½ m V² τ = −30 + 18
    assert s.shift == pytest.approx(s.shift_theory, rel=1e-8)  # どの経路でも同じ定数
    # 同じ形（定数ずれだけ）なので、停留点の位置（ε = 0）は K でも K' でも同じ
    mid = len(eps) // 2
    assert s.action_k - s.action_k[mid] == pytest.approx(s.action_kprime - s.action_kprime[mid])
    assert np.argmin(s.action_k) == np.argmin(s.action_kprime) == mid


def test_action_shift_vanishes_for_zero_boost():
    s = action_shift_under_boost(1.0, 0.0, 0.0, 0.0, 1.0, 2.0, np.linspace(-1, 1, 5))
    assert s.shift == pytest.approx(0.0, abs=1e-12)


@pytest.fixture(scope="module")
def frames():
    params = TwoBodyParams()
    return simulate_two_body_frames(params, np.array([0.35, 0.2]), duration=params.period)


def test_two_frames_give_the_same_orbit(frames):
    """K で独立に解いた軌道 = K' の軌道 + V t（運動方程式がガリレイ変換で形を変えない）。"""
    assert frames.frame_difference.max() < 1e-9


def test_momentum_and_center_of_mass_in_the_lab_frame(frames):
    mass = frames.params.total_mass
    assert frames.momentum_lab == pytest.approx(
        np.tile(mass * frames.frame_velocity, (frames.t.size, 1)), abs=1e-12
    )
    # 重心は等速直線運動: R(t) = V t
    assert frames.center_lab == pytest.approx(frames.t[:, None] * frames.frame_velocity, abs=1e-9)


def test_energy_is_conserved_in_each_frame_and_shifted_by_half_m_v_squared(frames):
    e0 = frames.energy_cm[0]
    assert np.max(np.abs(frames.energy_cm - e0)) / abs(e0) < 1e-5
    shift = 0.5 * frames.params.total_mass * float(frames.frame_velocity @ frames.frame_velocity)
    assert frames.energy_lab - frames.energy_cm == pytest.approx(shift, rel=1e-9)


def test_orbit_closes_after_one_kepler_period(frames):
    """ケプラーの第 3 法則の周期 T = 2π√2 で、相対位置が元に戻る。"""
    rel = frames.q_cm[:, 0:2] - frames.q_cm[:, 2:4]
    assert rel[-1] == pytest.approx(rel[0], abs=5e-3)
    separation = np.linalg.norm(rel, axis=1)
    assert separation.min() == pytest.approx(1.0, abs=1e-3)  # 近点
    assert separation.max() == pytest.approx(3.0, abs=1e-3)  # 遠点 = a(1+e) = 3
