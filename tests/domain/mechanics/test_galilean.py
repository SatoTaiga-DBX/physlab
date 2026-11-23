"""ガリレイ変換: 位置・速度の変換と、自由粒子の作用のずれ。"""

import numpy as np
import pytest

from physlab.domain.mechanics.galilean import (
    free_particle_action_shift,
    to_moving_frame_position,
    to_moving_frame_velocity,
)


def test_position_transform_1d_and_2d():
    t = np.array([0.0, 1.0, 2.0])
    x = np.array([5.0, 6.0, 7.0])
    assert to_moving_frame_position(x, t, 2.0) == pytest.approx([5.0, 4.0, 3.0])
    r = np.array([[0.0, 1.0], [1.0, 1.0], [2.0, 1.0]])
    shifted = to_moving_frame_position(r, t, np.array([1.0, 0.0]))
    assert shifted == pytest.approx(np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]))  # 一緒に動く


def test_velocity_transform_is_subtraction():
    assert to_moving_frame_velocity(np.array([3.0, 4.0]), np.array([1.0, 1.0])) == pytest.approx(
        [2.0, 3.0]
    )


def test_inverse_transform_recovers_original():
    t = np.linspace(0, 1, 5)
    x = np.sin(t)
    v_frame = 1.7
    assert to_moving_frame_position(x, t, v_frame) + v_frame * t == pytest.approx(x)


def test_action_shift_special_cases():
    assert free_particle_action_shift(2.0, 0.0, 1.0, 3.0) == 0.0  # 同じ系ならずれない
    assert free_particle_action_shift(2.0, 3.0, 0.0, 1.0) == pytest.approx(0.5 * 2.0 * 9.0 * 1.0)
