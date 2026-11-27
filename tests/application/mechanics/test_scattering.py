"""散乱の数える実験と、理論の断面積の照合。"""

import numpy as np
import pytest

from physlab.application.mechanics.deflection import deflection_in_chunks
from physlab.application.mechanics.scattering import (
    branch_cross_section,
    count_cross_section,
    curve_cross_section,
    deflection_curve,
    hard_sphere_angle_from_reflection,
    hard_sphere_lab_cross_section,
    hill,
    rainbow,
    rainbow_edges,
    solid_angles,
)

ENERGY = 1.5


def test_solid_angles_cover_the_sphere() -> None:
    assert solid_angles(np.linspace(0.0, np.pi, 7)).sum() == pytest.approx(4 * np.pi)


def test_hard_sphere_cross_section_is_isotropic() -> None:
    edges = np.linspace(0.05, np.pi, 19)
    counted = count_cross_section(
        lambda rho: hard_sphere_angle_from_reflection(rho, 1.0), 1.0, edges
    )
    assert np.max(np.abs(counted.measured - 0.25) / counted.error) < 4.0
    total = np.sum(counted.measured * solid_angles(edges))
    assert total == pytest.approx(np.pi, rel=0.01)  # 全断面積 π a²


@pytest.fixture(scope="module")
def curve():
    return deflection_curve(hill(), ENERGY, 3.5)


def test_hill_has_a_rainbow_angle(curve) -> None:
    rho_r, chi_r = rainbow(curve)
    assert 0.2 < rho_r < 0.7
    assert np.degrees(chi_r) == pytest.approx(36.2, abs=0.2)
    assert curve.chi.max() == pytest.approx(chi_r, abs=1e-5)


def test_counting_agrees_with_the_deflection_function(curve) -> None:
    _, chi_r = rainbow(curve)
    edges = np.radians(np.linspace(2.0, np.degrees(chi_r) + 3.0, 25))
    counted = count_cross_section(
        lambda rho: deflection_in_chunks(rho, hill(), ENERGY), 3.5, edges, seed=2
    )
    theory = curve_cross_section(curve, edges)
    used = theory > 0
    assert np.max(np.abs(counted.measured - theory)[used] / counted.error[used]) < 4.0
    assert counted.counts[~used].sum() == 0  # 虹の角より外には飛ばない


@pytest.mark.parametrize("energy", [1.05, 1.5, 30.0])
def test_counting_range_covers_the_rainbow_at_any_energy(energy: float) -> None:
    """高いエネルギーでは虹の角が 2° より小さくなる。区間が空にならないこと。"""
    fine = deflection_curve(hill(), energy, 3.5)
    _, chi_r = rainbow(fine)
    edges = rainbow_edges(chi_r)
    assert edges[0] < chi_r < edges[-1]
    counted = count_cross_section(
        lambda rho: deflection_in_chunks(rho, hill(), energy), 3.5, edges, samples=50_000
    )
    theory = curve_cross_section(fine, edges)
    used = theory > 0
    assert used.sum() > 10
    assert np.max(np.abs(counted.measured - theory)[used] / counted.error[used]) < 4.5


def test_counting_range_at_the_textbook_energy_is_two_degrees_to_rainbow_plus_three() -> None:
    edges = np.degrees(rainbow_edges(np.radians(36.25)))
    assert (edges[0], edges[-1]) == pytest.approx((2.0, 39.25))


def test_branch_formula_sums_two_impact_parameters(curve) -> None:
    """虹の角の手前では、2 つの ρ が同じ角に出る。角の平均は面積の足し合わせと合う。"""
    chi = np.radians(np.linspace(20.0, 30.0, 201))
    pointwise = branch_cross_section(curve, chi)
    edges = np.radians([20.0, 30.0])
    omega = 2 * np.pi * np.sin(chi)
    average = np.trapezoid(pointwise * omega, chi) / np.trapezoid(omega, chi)
    assert average == pytest.approx(curve_cross_section(curve, edges)[0], rel=2e-3)


@pytest.mark.parametrize(("m1", "m2"), [(1.0, 5.0), (1.0, 1.0), (5.0, 1.0)])
def test_lab_frame_of_hard_spheres_follows_the_decay_distribution(m1: float, m2: float) -> None:
    result = hard_sphere_lab_cross_section(m1, m2)
    used = result.expected > 0
    z = np.abs(result.counted.measured - result.expected)[used] / result.counted.error[used]
    assert np.max(z) < 4.0
    assert result.counted.counts.sum() == result.counted.samples  # すべて上限の角の内側
