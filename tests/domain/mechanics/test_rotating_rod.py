"""回転する棒の上の質点: 解析解が運動方程式を満たし、E = pṙ − L は一定、T は一定でない。"""

import numpy as np
import pytest

from physlab.domain.mechanics.rotating_rod import (
    energy_function,
    kinetic_energy,
    motor_power,
    radius,
)


def test_solution_satisfies_the_equation_of_motion() -> None:
    omega, r0 = 1.7, 0.2
    t = np.linspace(0.0, 1.5, 3001)
    r, _ = radius(t, r0, omega)
    dt = t[1] - t[0]
    r_ddot = (r[2:] - 2 * r[1:-1] + r[:-2]) / dt**2
    assert r_ddot == pytest.approx(omega**2 * r[1:-1], rel=1e-5)


def test_h_is_constant_but_kinetic_energy_is_not() -> None:
    omega, r0, mass = 1.5, 0.1, 2.0
    t = np.linspace(0.0, 2.0, 50)
    r, r_dot = radius(t, r0, omega)
    h = energy_function(r, r_dot, omega, mass)
    assert h == pytest.approx(-0.5 * mass * omega**2 * r0**2, rel=1e-9)
    kinetic = kinetic_energy(r, r_dot, omega, mass)
    assert np.all(np.diff(kinetic) > 0)
    # T − h = m ω² r²（遠心項のぶんだけずれる）
    assert kinetic - h == pytest.approx(mass * omega**2 * r**2, rel=1e-12)


def test_motor_power_is_the_rate_of_change_of_kinetic_energy() -> None:
    omega, r0 = 1.3, 0.3
    t, s = 0.9, 1e-6
    k = [float(kinetic_energy(*radius(t + d, r0, omega), omega)) for d in (-s, s)]
    power = float(motor_power(*radius(t, r0, omega), omega))
    assert (k[1] - k[0]) / (2 * s) == pytest.approx(power, rel=1e-6)
