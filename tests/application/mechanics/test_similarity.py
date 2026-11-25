"""力学的相似: 周期の冪、相似変換した軌道の重なり、ビリアル定理、ケプラーの第 3 法則。"""

import pytest

from physlab.application.mechanics.similarity import (
    collapse,
    kepler_periods,
    kepler_virial,
    sweep_periods,
    virial,
)


@pytest.fixture(scope="module")
def sweep():
    return sweep_periods()


@pytest.mark.parametrize(("k", "rel"), [(1.0, 2e-3), (2.0, 1e-6), (4.0, 1e-6), (6.0, 1e-6)])
def test_measured_periods_match_the_exact_formula(sweep, k: float, rel: float) -> None:
    # k = 1 は原点で力が不連続なので、数値解の精度が落ちる
    assert sweep.measured[k] == pytest.approx(sweep.exact[k], rel=rel)
    assert sweep.slopes[k] == pytest.approx(1 - k / 2, abs=1e-3)


def test_kepler_periods_scale_as_the_three_halves_power() -> None:
    kepler = kepler_periods()
    assert kepler.measured == pytest.approx(kepler.exact, rel=1e-5)
    assert kepler.slope == pytest.approx(1.5, abs=1e-5)


@pytest.mark.parametrize("k", [1.0, 4.0, 6.0])
def test_similar_motions_collapse_onto_one_curve(k: float) -> None:
    result = collapse(k)
    assert result.max_difference < 1e-9


def test_virial_theorem_over_whole_periods() -> None:
    result = virial()
    for k, ratio in result.over_whole_periods.items():
        assert ratio == pytest.approx(k, rel=1e-5), k


def test_kepler_virial_ratio_is_minus_one() -> None:
    assert kepler_virial() == pytest.approx(-1.0, abs=1e-5)
