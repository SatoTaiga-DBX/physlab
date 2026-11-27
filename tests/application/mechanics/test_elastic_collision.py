"""2 粒子の運動方程式を解いた衝突が、保存則の予言する円に乗ること。"""

import numpy as np
import pytest

from physlab.application.mechanics.elastic_collision import (
    DEMO_POTENTIALS,
    circle_experiment,
    outcome,
    simulate_collisions,
)
from physlab.domain.mechanics.collision import (
    energy_transfer_fraction,
    lab_angles,
    momentum_circle,
)

IMPACT = np.array([0.0, 0.2, 0.45, 0.7, 0.9])


@pytest.fixture(scope="module")
def results():
    return circle_experiment(1.0, 2.0, IMPACT)


def test_final_momenta_lie_on_the_same_circle_for_every_force(results) -> None:
    center, radius = momentum_circle(1.0, 2.0, 1.0)
    assert set(results) == set(DEMO_POTENTIALS)
    for name, result in results.items():
        distance = np.hypot(result.p2[:, 0] - center, result.p2[:, 1])
        assert distance == pytest.approx(radius, abs=1e-5), name
        assert result.p1 + result.p2 == pytest.approx(np.tile([1.0, 0.0], (5, 1)), abs=1e-10)


def test_lab_angles_and_energy_transfer_follow_from_chi(results) -> None:
    for name, result in results.items():
        theta1, theta2 = lab_angles(result.chi, 1.0, 2.0)
        assert result.theta1 == pytest.approx(theta1, abs=1e-4), name
        hit = result.chi > 1e-3
        assert result.theta2[hit] == pytest.approx(theta2[hit], abs=1e-4), name
        assert result.transfer == pytest.approx(
            energy_transfer_fraction(result.chi, 1.0, 2.0), abs=1e-4
        )


def test_the_force_decides_only_where_on_the_circle(results) -> None:
    """同じ ρ でも、力の形が違えば χ は違う。"""
    chi = np.array([result.chi[2] for result in results.values()])
    assert chi.max() - chi.min() > 0.3


def test_head_on_collision_reflects_in_the_center_of_mass_frame(results) -> None:
    for result in results.values():
        assert result.chi[0] == pytest.approx(np.pi, abs=1e-9)


def test_equal_masses_leave_at_right_angles() -> None:
    run = simulate_collisions(1.0, 1.0, DEMO_POTENTIALS["gaussian"], [0.3, 0.8])
    result = outcome(run)
    assert result.theta1 + result.theta2 == pytest.approx(np.pi / 2, abs=1e-9)


def test_energy_is_conserved_during_the_collision() -> None:
    run = simulate_collisions(1.0, 2.0, DEMO_POTENTIALS["steep-core"], IMPACT)
    assert np.max(np.abs(run.energy - run.energy[0])) < 1e-4
