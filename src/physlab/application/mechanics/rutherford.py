"""ラザフォード散乱の実験: クーロン力で打ち込んで数える。電荷を広げた場合と比べる。

設計書: docs/design/007-collisions.md ／ 教科書: 第 5 章 §ラザフォードの公式

単位は、斥力の正面衝突の最接近距離 ``d = α/E`` が 1 になるように取る
（``mass = speed = 1``、``α = 0.5``）。偏向角は公式を使わず、軌道の式の求積で求める。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from physlab.application.mechanics.deflection import deflection, deflection_in_chunks
from physlab.application.mechanics.scattering import (
    CrossSection,
    DeflectionCurve,
    count_cross_section,
)
from physlab.domain.mechanics.collision import final_velocities
from physlab.domain.mechanics.decay import reduced_mass
from physlab.domain.mechanics.potentials import UniformBall, coulomb
from physlab.domain.mechanics.scattering import (
    energy_transfer_cross_section,
    max_energy_transfer,
    rutherford_impact_parameter,
)
from physlab.numerics.integrators import Array
from physlab.numerics.quadrature import integrate_smooth

ALPHA = 0.5
MASS = 1.0
SPEED = 1.0
ENERGY = 0.5 * MASS * SPEED**2
HEAD_ON = ALPHA / ENERGY  # d = 1


def point_charge_cross_section(
    min_angle: float = np.radians(20.0),
    n_bins: int = 16,
    samples: int = 200_000,
    seed: int = 1,
) -> CrossSection:
    """``min_angle`` より大きく曲がる粒子を打ち込んで数える。"""
    # ρ_max は min_angle をわずかに下回る角まで届くように取る
    rho_max = float(rutherford_impact_parameter(0.9 * min_angle, HEAD_ON))
    edges = np.linspace(min_angle, np.pi, n_bins + 1)
    return count_cross_section(
        lambda rho: deflection_in_chunks(rho, coulomb(ALPHA), ENERGY),
        rho_max,
        edges,
        samples=samples,
        seed=seed,
    )


def ball_cross_section(
    radius_over_d: float,
    min_angle: float = np.radians(20.0),
    n_bins: int = 16,
    samples: int = 200_000,
    seed: int = 1,
) -> CrossSection:
    """電荷が半径 ``R = radius_over_d × d`` の球に一様に広がっている場合。"""
    potential = UniformBall(ALPHA, radius_over_d * HEAD_ON)
    rho_max = max(
        float(rutherford_impact_parameter(0.9 * min_angle, HEAD_ON)), 1.01 * potential.radius
    )
    edges = np.linspace(min_angle, np.pi, n_bins + 1)
    return count_cross_section(
        lambda rho: deflection_in_chunks(rho, potential, ENERGY),
        rho_max,
        edges,
        samples=samples,
        seed=seed,
    )


def ball_curve(radius_over_d: float, rho_max: float = 4.0, n: int = 1500) -> DeflectionCurve:
    """帯電球の偏向角の曲線。表面をかすめる付近を含めて等間隔に取る。"""
    potential = UniformBall(ALPHA, radius_over_d * HEAD_ON)
    rho = np.linspace(rho_max / n, rho_max, n)
    return DeflectionCurve(rho, deflection(rho, potential, ENERGY))


def max_ball_angle(radius_over_d: float, n: int = 2000) -> float:
    """帯電球（``R > 1.5 d``、中心を通り抜けられる）の散乱角の最大値。"""
    radius = radius_over_d * HEAD_ON
    rho = np.linspace(radius / n, 1.2 * radius, n)
    return float(np.max(np.abs(deflection(rho, UniformBall(ALPHA, radius), ENERGY))))


@dataclass(frozen=True)
class EnergyTransfer:
    """静止した標的に渡るエネルギー ε の分布 ``dσ/dε``。"""

    m1: float
    m2: float
    maximum: float  # 渡せる最大のエネルギー
    edges: Array
    measured: Array
    error: Array
    expected: Array  # 区間の平均

    @property
    def centers(self) -> Array:
        return np.sqrt(self.edges[:-1] * self.edges[1:])


def energy_transfer(
    m1: float = 1.0,
    m2: float = 2.0,
    min_fraction: float = 0.05,
    n_bins: int = 14,
    samples: int = 200_000,
    seed: int = 1,
) -> EnergyTransfer:
    """換算質量の運動として打ち込み、実験室系の標的の速度から ε を数える。"""
    mu = reduced_mass(m1, m2)
    energy = 0.5 * mu * SPEED**2
    d = ALPHA / energy
    maximum = max_energy_transfer(m1, m2, SPEED)
    # ε = ε_max sin²(χ/2) が下限になる角より少し外側まで打ち込む
    lowest = 2.0 * np.arcsin(np.sqrt(min_fraction))
    rho_max = float(rutherford_impact_parameter(0.9 * lowest, d))
    rng = np.random.default_rng(seed)
    rho = np.maximum(rho_max * np.sqrt(rng.uniform(0.0, 1.0, samples)), 1e-9 * rho_max)
    chi = deflection_in_chunks(rho, coulomb(ALPHA), energy)
    direction = np.column_stack([np.cos(chi), np.sin(chi)])
    _, u2 = final_velocities(m1, m2, [SPEED, 0.0], [0.0, 0.0], direction)
    epsilon = 0.5 * m2 * np.sum(u2**2, axis=1)

    edges = np.geomspace(min_fraction * maximum, maximum, n_bins + 1)
    edges[-1] *= 1.0 + 1e-9
    counts, _ = np.histogram(epsilon, bins=edges)
    area = np.pi * rho_max**2 / samples
    width = np.diff(edges)
    expected = integrate_smooth(
        lambda e: energy_transfer_cross_section(e, ALPHA, m2, SPEED), edges[:-1], edges[1:]
    )
    expected = expected / width
    return EnergyTransfer(
        m1=m1,
        m2=m2,
        maximum=maximum,
        edges=edges,
        measured=area * counts / width,
        error=area * np.sqrt(np.maximum(counts, 1)) / width,
        expected=expected,
    )
