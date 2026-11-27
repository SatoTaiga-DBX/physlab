"""偏向角の求積を、解析解と、運動方程式の時間積分の両方と照合する。"""

import numpy as np
import pytest

from physlab.application.mechanics.deflection import (
    closest_distance,
    deflection,
    incoming_state,
    orbit_deflection,
    trace_orbits,
    turning_u,
)
from physlab.application.mechanics.scattering import hill
from physlab.domain.mechanics.potentials import PowerLaw, UniformBall, coulomb
from physlab.domain.mechanics.scattering import (
    coulomb_closest_approach,
    inverse_square_deflection,
    rutherford_angle,
)

RHO = np.array([0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 20.0])


@pytest.mark.parametrize("sign", [1.0, -1.0], ids=["repulsive", "attractive"])
def test_coulomb_quadrature_matches_rutherford(sign: float) -> None:
    theta = deflection(RHO, coulomb(sign * 0.5), 0.5)
    assert theta == pytest.approx(sign * rutherford_angle(RHO, 1.0), abs=1e-8)


def test_inverse_square_quadrature_matches_the_exact_solution() -> None:
    theta = deflection(RHO, PowerLaw(0.5, 2.0), 0.5)
    assert theta == pytest.approx(inverse_square_deflection(RHO, 0.5, 0.5), abs=1e-8)


def test_turning_point_satisfies_the_energy_equation() -> None:
    u = turning_u(RHO, coulomb(0.5), 0.5)
    chi = rutherford_angle(RHO, 1.0)
    assert 1.0 / u == pytest.approx(coulomb_closest_approach(chi, 1.0), rel=1e-10)


def test_rejects_zero_impact_parameter() -> None:
    with pytest.raises(ValueError):
        turning_u(0.0, coulomb(0.5), 0.5)


def test_incoming_state_has_the_asymptotic_energy_and_angular_momentum() -> None:
    rho = np.array([0.5, 2.0])
    q, v = incoming_state(coulomb(0.5), 1.0, 1.0, rho, 12.0)
    x, y = q[:2], q[2:]
    vx, vy = v[:2], v[2:]
    r = np.hypot(x, y)
    assert r == pytest.approx(12.0)
    assert 0.5 * (vx**2 + vy**2) + 0.5 / r == pytest.approx(0.5)
    assert x * vy - y * vx == pytest.approx(-rho)


@pytest.mark.parametrize(
    ("potential", "energy"),
    [(hill(), 1.5), (UniformBall(0.5, 1.4), 0.5), (UniformBall(0.5, 3.0), 0.5)],
    ids=["gaussian", "ball-1.4", "ball-3"],
)
def test_time_integration_agrees_with_the_quadrature(potential, energy: float) -> None:
    rho = np.array([0.2, 0.7, 1.2, 2.5])
    orbits = trace_orbits(potential, 1.0, np.sqrt(2 * energy), rho, start_radius=10.0)
    assert orbit_deflection(orbits) == pytest.approx(deflection(rho, potential, energy), abs=1e-5)


@pytest.mark.parametrize("sign", [1.0, -1.0], ids=["repulsive", "attractive"])
def test_time_integration_of_coulomb_orbits(sign: float) -> None:
    rho = np.array([0.5, 1.0, 2.0])
    orbits = trace_orbits(coulomb(sign * 0.5), 1.0, 1.0, rho, start_radius=15.0)
    assert orbit_deflection(orbits) == pytest.approx(sign * rutherford_angle(rho, 1.0), abs=1e-5)
    if sign > 0:
        expected = coulomb_closest_approach(rutherford_angle(rho, 1.0), 1.0)
        assert closest_distance(orbits) == pytest.approx(expected, abs=1e-4)


def test_negative_impact_parameter_is_the_mirror_image() -> None:
    orbits = trace_orbits(hill(), 1.0, 2.0, [0.5, -0.5], start_radius=6.0)
    assert orbits.y[:, 1] == pytest.approx(-orbits.y[:, 0])
    theta = orbit_deflection(orbits)
    assert theta[1] == pytest.approx(-theta[0])
