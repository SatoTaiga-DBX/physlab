"""崩壊の乱数実験: 保存則と、理論の分布との照合。"""

import numpy as np
import pytest

from physlab.application.mechanics.decay import (
    angle_histogram,
    energy_error,
    energy_histogram,
    momentum_error,
    planar_burst,
    simulate_decays,
)
from physlab.domain.mechanics.decay import DecayParams

PARAMS = DecayParams()  # v10 = 6


@pytest.fixture(scope="module", params=[3.0, 6.0, 12.0], ids=["k=0.5", "k=1", "k=2"])
def sample(request):
    return simulate_decays(PARAMS, request.param, n=200_000, seed=3)


def test_every_decay_conserves_momentum_and_energy(sample) -> None:
    assert momentum_error(sample) < 1e-12
    assert energy_error(sample) < 1e-10


def test_angle_distribution_agrees_with_the_formula(sample) -> None:
    h = angle_histogram(sample)
    assert h.measured.sum() == pytest.approx(1.0)
    assert h.expected.sum() == pytest.approx(1.0, abs=1e-8)
    z = (h.measured - h.expected) / h.error
    assert np.max(np.abs(z)) < 4.0


def test_lab_energy_is_uniformly_distributed(sample) -> None:
    h = energy_histogram(sample)
    assert h.measured.sum() == pytest.approx(1.0)
    z = (h.measured - h.expected) / h.error
    assert np.max(np.abs(z)) < 4.0


def test_planar_burst_lies_on_the_rest_frame_circle() -> None:
    velocity = planar_burst(PARAMS, 24)
    assert np.linalg.norm(velocity, axis=1) == pytest.approx(6.0)
    assert velocity.mean(axis=0) == pytest.approx([0.0, 0.0], abs=1e-12)
