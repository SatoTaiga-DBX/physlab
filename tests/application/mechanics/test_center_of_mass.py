"""慣性中心: 投げた系の重心は放物線を描き、内部エネルギーと重心まわりの角運動量は保存する。"""

import numpy as np
import pytest

from physlab.application.mechanics.center_of_mass import throw_cluster


@pytest.fixture(scope="module")
def throw():
    return throw_cluster()


def test_center_of_mass_follows_a_parabola_despite_spinning_and_vibrating(throw) -> None:
    assert throw.center == pytest.approx(throw.center_exact, abs=1e-10)
    # 内部は実際に回って振動している
    assert np.ptp(throw.spring) > 1.0
    assert throw.center[0] == pytest.approx(0.0, abs=1e-12)
    assert abs(throw.center[-1, 1]) < 0.05  # 投げた高さに戻る


def test_internal_energy_and_spin_angular_momentum_are_conserved(throw) -> None:
    internal = throw.internal_energy
    assert np.max(np.abs(internal / internal[0] - 1)) < 1e-5
    spin = throw.angular_relative
    assert np.max(np.abs(spin - spin[0])) < 1e-10


def test_energy_splits_into_center_and_internal_parts(throw) -> None:
    total = throw.total_energy
    assert np.max(np.abs(total / total[0] - 1)) < 1e-5
    # 重心の運動エネルギー + 重力ポテンシャル は単独でも保存する
    center_part = throw.kinetic_center + throw.gravity
    assert center_part == pytest.approx(center_part[0], rel=1e-9)


def test_angular_momentum_splits_and_the_center_part_changes_by_gravity_torque(throw) -> None:
    assert throw.angular_total == pytest.approx(
        throw.angular_center + throw.angular_relative, abs=1e-9
    )
    assert abs(throw.angular_total[-1] - throw.angular_total[0]) > 100
