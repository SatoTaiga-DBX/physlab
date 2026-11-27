"""微小角散乱の実験: 直線の経路に沿った力積の近似と、厳密な偏向角を比べる。

設計書: docs/design/007-collisions.md ／ 教科書: 第 5 章 §微小角の散乱

力積の近似は ``θ = −(2ρ / m v²) ∫_ρ^∞ U'(r) dr / √(r² − ρ²)``。
``r = ρ / cos u`` と置くと端の特異性が消えて ``∫0^{π/2} U'(ρ/cos u) / cos u du`` になる。
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from physlab.application.mechanics.deflection import deflection, trace_orbits
from physlab.domain.mechanics.potentials import CentralPotential, Gaussian, PowerLaw
from physlab.numerics.integrators import Array
from physlab.numerics.quadrature import integrate_smooth

MASS = 1.0
SPEED = 1.0
ENERGY = 0.5 * MASS * SPEED**2


def impulse_angle(
    rho: ArrayLike,
    potential: CentralPotential,
    mass: float = MASS,
    speed: float = SPEED,
    n: int = 120,
) -> Array:
    """直線 ``y = ρ`` に沿って受ける横向きの力積から求めた偏向角。"""
    rho = np.atleast_1d(np.asarray(rho, dtype=float))
    column = rho[:, None]

    def integrand(u: Array) -> Array:
        return potential.derivative(column / np.cos(u)) / np.cos(u)

    split = np.full_like(rho, 0.25 * np.pi)
    if potential.breakpoint_radius is not None:
        ratio = rho / potential.breakpoint_radius
        split = np.where(ratio < 1.0, np.arccos(np.minimum(ratio, 1.0)), split)
    total = integrate_smooth(integrand, np.zeros_like(rho), split, n=n)
    total += integrate_smooth(integrand, split, np.full_like(rho, 0.5 * np.pi), n=n)
    return -2.0 * rho / (mass * speed**2) * total


# 強さ s を変えられる 4 種類の力
FAMILIES: dict[str, Callable[[float], CentralPotential]] = {
    "coulomb": lambda s: PowerLaw(alpha=s, n=1.0),
    "inverse-square": lambda s: PowerLaw(alpha=s, n=2.0),
    "inverse-fourth": lambda s: PowerLaw(alpha=s, n=4.0),
    "gaussian": lambda s: Gaussian(height=s, width=1.0),
}


@dataclass(frozen=True)
class ErrorSweep:
    """同じ ρ で力の強さを変えたときの、力積の近似の誤差。"""

    strength: Array
    exact: Array
    approx: Array

    @property
    def relative_error(self) -> Array:
        return np.abs(self.approx / self.exact - 1.0)

    @property
    def slope(self) -> float:
        """誤差が偏向角の何乗で減るか（両対数の傾き。角の小さい側の半分で測る）。"""
        half = self.exact.size // 2
        x = np.log(self.exact[:half])
        y = np.log(self.relative_error[:half])
        return float(np.polyfit(x, y, 1)[0])


def error_sweep(family: str, rho: float = 1.0, strengths: Array | None = None) -> ErrorSweep:
    strengths = np.geomspace(3e-3, 0.3, 13) if strengths is None else strengths
    make = FAMILIES[family]
    exact = np.array([deflection(rho, make(s), ENERGY)[0] for s in strengths])
    approx = np.array([impulse_angle(rho, make(s))[0] for s in strengths])
    return ErrorSweep(strengths, exact, approx)


@dataclass(frozen=True)
class CrossSectionTail:
    """べき型の力の、小さい角での ``dσ/dΩ``。"""

    n: float
    theta: Array
    cross_section: Array
    slope: float  # 両対数の傾き。理論は −2 − 2/n


def small_angle_tail(n: float, alpha: float = 0.5) -> CrossSectionTail:
    """遠くを通る ρ の厳密な偏向角から ``dσ/dΩ = (ρ / sin θ) |dρ/dθ|`` を作る。"""
    potential = PowerLaw(alpha=alpha, n=n)
    # 力積の近似で θ が 0.002〜0.05 になる ρ の範囲を選ぶ
    unit = float(impulse_angle(1.0, potential)[0])
    rho = np.geomspace((unit / 0.05) ** (1.0 / n), (unit / 0.002) ** (1.0 / n), 400)
    theta = deflection(rho, potential, ENERGY)
    cross_section = rho / np.sin(theta) / np.abs(np.gradient(theta, rho))
    inner = slice(2, -2)  # 端の数値微分は精度が落ちる
    slope = float(np.polyfit(np.log(theta[inner]), np.log(cross_section[inner]), 1)[0])
    return CrossSectionTail(n, theta[inner], cross_section[inner], slope)


@dataclass(frozen=True)
class ImpulseDemo:
    """ほとんど曲がらない 1 本の軌道と、それを直線とみなした近似。"""

    rho: float
    potential: CentralPotential
    t: Array
    x: Array
    y: Array
    force_y: Array  # 実際の軌道に沿った横向きの力
    straight_x: Array
    straight_force_y: Array  # 直線 y = ρ に沿った横向きの力
    momentum_y: Array  # 実際の横向きの運動量の時間変化
    approx_momentum_y: float  # 力積の近似の最終値
    exact_angle: float
    approx_angle: float


def impulse_demo(energy_ratio: float = 10.0, rho: float = 0.7, start: float = 5.0) -> ImpulseDemo:
    """高さ 1 のガウス型の丘を、エネルギー ``energy_ratio`` で通り抜ける。"""
    potential = Gaussian(height=1.0, width=1.0)
    speed = float(np.sqrt(2.0 * energy_ratio / MASS))
    orbits = trace_orbits(
        potential, MASS, speed, [rho], start_radius=start, duration=2 * start / speed, dt=5e-4
    )
    x, y = orbits.x[:, 0], orbits.y[:, 0]
    r = np.hypot(x, y)
    force_y = -potential.derivative(r) * y / r
    straight = np.linspace(-start, start, 801)
    sr = np.hypot(straight, rho)
    straight_force = -potential.derivative(sr) * rho / sr
    approx = float(impulse_angle(rho, potential, MASS, speed)[0])
    exact = float(deflection(rho, potential, 0.5 * MASS * speed**2)[0])
    return ImpulseDemo(
        rho=rho,
        potential=potential,
        t=orbits.t,
        x=x,
        y=y,
        force_y=force_y,
        straight_x=straight,
        straight_force_y=straight_force,
        momentum_y=MASS * orbits.vy[:, 0],
        approx_momentum_y=approx * MASS * speed,
        exact_angle=exact,
        approx_angle=approx,
    )
