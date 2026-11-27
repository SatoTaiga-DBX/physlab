"""散乱の公式どうしの整合性。"""

import numpy as np
import pytest

from physlab.domain.mechanics.potentials import coulomb
from physlab.domain.mechanics.scattering import (
    ball_graze_angle,
    ball_small_angle,
    coulomb_closest_approach,
    energy_transfer_cross_section,
    hard_sphere_angle,
    hard_sphere_reflect,
    head_on_distance,
    inverse_square_deflection,
    max_energy_transfer,
    power_law_small_angle,
    power_law_small_angle_cross_section,
    radial_function,
    rutherford_angle,
    rutherford_cross_section,
    rutherford_impact_parameter,
)

CHI = np.linspace(0.2, 3.0, 8)


def test_head_on_distance_is_where_the_particle_turns_back() -> None:
    d = head_on_distance(0.5, mass=1.0, speed=1.0)
    assert d == pytest.approx(1.0)
    f = radial_function(np.array(1.0 / d), np.array(0.0), coulomb(0.5), 0.5)
    assert float(f) == pytest.approx(0.0, abs=1e-15)


def test_rutherford_angle_and_impact_parameter_are_inverse() -> None:
    assert rutherford_angle(rutherford_impact_parameter(CHI, 2.0), 2.0) == pytest.approx(CHI)


def test_cross_section_is_rho_over_sin_chi_times_drho_dchi() -> None:
    h = 1e-6
    rho = rutherford_impact_parameter(CHI, 1.0)
    slope = (
        rutherford_impact_parameter(CHI + h, 1.0) - rutherford_impact_parameter(CHI - h, 1.0)
    ) / (2 * h)
    assert rutherford_cross_section(CHI, 1.0) == pytest.approx(
        rho / np.sin(CHI) * np.abs(slope), rel=1e-8
    )


def test_closest_approach_is_the_turning_point() -> None:
    rho = rutherford_impact_parameter(CHI, 1.0)
    r = coulomb_closest_approach(CHI, 1.0)
    assert radial_function(1.0 / r, rho, coulomb(0.5), 0.5) == pytest.approx(0.0, abs=1e-12)


def test_energy_transfer_cross_section_follows_from_the_angle_formula() -> None:
    """dσ = 2πρ|dρ| を ε = ε_max sin²(χ/2) で書き直す。"""
    m1, m2, speed, alpha = 1.0, 2.0, 1.0, 0.5
    mu = m1 * m2 / (m1 + m2)
    d = head_on_distance(alpha, mu, speed)
    top = max_energy_transfer(m1, m2, speed)
    lo, hi = 0.2 * top, 0.3 * top
    chi = 2 * np.arcsin(np.sqrt(np.array([lo, hi]) / top))
    area = np.pi * np.diff(rutherford_impact_parameter(chi, d) ** 2)[0]
    e = np.linspace(lo, hi, 20_001)
    integral = np.trapezoid(energy_transfer_cross_section(e, alpha, m2, speed), e)
    assert integral == pytest.approx(-area, rel=1e-8)


def test_hard_sphere_reflection_matches_the_angle_formula() -> None:
    rho = np.linspace(0.0, 0.99, 12)
    vx, vy = hard_sphere_reflect(rho, 1.0)
    assert np.hypot(vx, vy) == pytest.approx(1.0)
    assert np.arctan2(vy, vx) == pytest.approx(hard_sphere_angle(rho, 1.0), abs=1e-12)
    vx, vy = hard_sphere_reflect(np.array([1.5]), 1.0)  # 当たらないので直進
    assert (vx[0], vy[0]) == (1.0, 0.0)


def test_inverse_square_deflection_limits() -> None:
    assert inverse_square_deflection(0.0, 0.5, 0.5) == pytest.approx(np.pi)
    assert inverse_square_deflection(1e4, 0.5, 0.5) == pytest.approx(0.0, abs=1e-7)


@pytest.mark.parametrize("n", [1.0, 2.0, 4.0])
def test_small_angle_formulas(n: float) -> None:
    rho = np.array([50.0, 80.0])
    theta = power_law_small_angle(rho, 0.5, n, mass=1.0, speed=1.0)
    coefficient = theta[0] * rho[0] ** n
    assert theta[1] == pytest.approx(coefficient / rho[1] ** n)
    # dσ/dΩ ≈ (ρ/θ)|dρ/dθ|
    h = 1e-6 * theta[0]
    rho_of = lambda t: (coefficient / t) ** (1 / n)  # noqa: E731
    slope = (rho_of(theta[0] + h) - rho_of(theta[0] - h)) / (2 * h)
    assert power_law_small_angle_cross_section(theta[0], coefficient, n) == pytest.approx(
        rho[0] / theta[0] * abs(slope), rel=1e-7
    )


def test_coulomb_small_angle_agrees_with_rutherford() -> None:
    rho = 1000.0
    assert power_law_small_angle(rho, 0.5, 1.0, 1.0, 1.0) == pytest.approx(
        rutherford_angle(rho, 1.0), rel=1e-6
    )
    assert power_law_small_angle(rho, 0.5, 2.0, 1.0, 1.0) == pytest.approx(
        inverse_square_deflection(rho, 0.5, 0.5), rel=1e-5
    )


def test_ball_small_angle_is_the_point_charge_outside_and_weaker_inside() -> None:
    assert ball_small_angle(5.0, 1.0, 2.0) == pytest.approx(1.0 / 5.0)
    assert ball_small_angle(1.0, 1.0, 2.0) == pytest.approx(1.0 - 0.75**1.5)


def test_ball_graze_angle() -> None:
    assert ball_graze_angle(0.8) == pytest.approx(np.pi)
    chi = ball_graze_angle(1.4)
    assert coulomb_closest_approach(chi, 1.0) == pytest.approx(1.4)
