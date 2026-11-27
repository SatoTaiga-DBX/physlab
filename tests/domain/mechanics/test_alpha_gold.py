"""α 粒子と金の実験の大きさ。"""

import numpy as np
import pytest

from physlab.domain.mechanics.alpha_gold import AlphaOnGold
from physlab.domain.mechanics.collision import max_deflection


def test_head_on_distance_is_much_larger_than_the_nucleus() -> None:
    setup = AlphaOnGold()
    assert setup.head_on_fm == pytest.approx(41.37, abs=0.01)
    assert setup.nuclear_radius_fm == pytest.approx(6.98, abs=0.01)
    assert setup.head_on_fm > 5 * setup.nuclear_radius_fm


def test_an_electron_barely_deflects_an_alpha_particle() -> None:
    assert AlphaOnGold().max_angle_from_electron == pytest.approx(
        max_deflection(7294.3, 1.0), rel=1e-12
    )


def test_fraction_beyond_ninety_degrees() -> None:
    setup = AlphaOnGold()
    sigma = np.pi * (0.5 * setup.head_on_fm * 1e-15) ** 2
    assert setup.cross_section_beyond_m2(np.pi / 2) == pytest.approx(sigma)
    assert setup.fraction_beyond(np.pi / 2) == pytest.approx(5.9e28 * 4e-7 * sigma)
