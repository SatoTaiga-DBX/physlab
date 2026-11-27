"""散乱の実験: 偏向角の曲線を求め、乱数で打ち込んで角度ごとに数える。

設計書: docs/design/007-collisions.md ／ 教科書: 第 5 章 §散乱

数える実験では、入射方向に垂直な半径 ``ρ_max`` の円板内に一様に粒子を打ち込む
（``ρ = ρ_max √ξ``、ξ は 0〜1 の一様乱数）。角度の区間に入った個数 ``ΔN`` から
``dσ/dΩ ≈ (π ρ_max² ΔN / N) / ΔΩ``。
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from physlab.application.mechanics.decay import angle_bin_probability
from physlab.application.mechanics.deflection import deflection
from physlab.domain.mechanics.collision import final_velocities
from physlab.domain.mechanics.decay import max_lab_angle
from physlab.domain.mechanics.potentials import CentralPotential, Gaussian
from physlab.domain.mechanics.scattering import hard_sphere_reflect
from physlab.numerics.integrators import Array


def hill(width: float = 1.0) -> Gaussian:
    """実験に使うガウス型の丘（高さ 1）。区間を分ける半径を中心近くに置く。"""
    return Gaussian(height=1.0, width=width, breakpoint_radius=width / 3.0)


@dataclass(frozen=True)
class DeflectionCurve:
    rho: Array
    theta: Array  # 符号つきの偏向角 Θ

    @property
    def chi(self) -> Array:
        return np.abs(self.theta)


@dataclass(frozen=True)
class CrossSection:
    """角度の区間ごとの ``dσ/dΩ``。"""

    edges: Array  # rad
    measured: Array
    error: Array  # 統計誤差（1σ）
    counts: Array
    samples: int

    @property
    def centers(self) -> Array:
        return 0.5 * (self.edges[:-1] + self.edges[1:])


def solid_angles(edges: Array) -> Array:
    """角度の区間がもつ立体角 ``2π (cos χ1 − cos χ2)``。"""
    return 2.0 * np.pi * (np.cos(edges[:-1]) - np.cos(edges[1:]))


def deflection_curve(
    potential: CentralPotential, energy: float, rho_max: float, n: int = 4000
) -> DeflectionCurve:
    rho = np.linspace(rho_max / n, rho_max, n)
    return DeflectionCurve(rho, deflection(rho, potential, energy))


def rainbow(curve: DeflectionCurve) -> tuple[float, float]:
    """偏向角が極大になる ``(ρ_r, χ_r)``（放物線で補間）。"""
    chi = curve.chi
    i = int(np.argmax(chi))
    if i in (0, chi.size - 1):
        raise ValueError("偏向角に内側の極大が無い")
    a, b, c = np.polyfit(curve.rho[i - 1 : i + 2], chi[i - 1 : i + 2], 2)
    peak = -b / (2.0 * a)
    return float(peak), float(a * peak**2 + b * peak + c)


def rainbow_edges(rainbow_angle: float, n_bins: int = 24) -> Array:
    """虹の角の手前から少し外側までを数える区間（rad）。

    ふつうは 2° から虹の角 + 3° まで。エネルギーが高くて虹の角が小さいときは、
    虹の角に比例して狭める（そうしないと、区間が虹の角の外だけになる）。
    """
    low = min(np.radians(2.0), 0.1 * rainbow_angle)
    high = rainbow_angle + min(np.radians(3.0), 0.1 * rainbow_angle)
    return np.linspace(low, high, n_bins + 1)


def count_cross_section(
    angle_of: Callable[[Array], Array],
    rho_max: float,
    edges: Array,
    samples: int = 200_000,
    seed: int = 1,
) -> CrossSection:
    """``angle_of(ρ)`` で散乱される粒子を乱数で打ち込み、角度ごとに数える。"""
    rng = np.random.default_rng(seed)
    rho = rho_max * np.sqrt(rng.uniform(0.0, 1.0, samples))
    rho = np.maximum(rho, 1e-9 * rho_max)
    chi = np.abs(angle_of(rho))
    counts, _ = np.histogram(chi, bins=edges)
    area_per_particle = np.pi * rho_max**2 / samples
    omega = solid_angles(edges)
    return CrossSection(
        edges=edges,
        measured=area_per_particle * counts / omega,
        error=area_per_particle * np.sqrt(np.maximum(counts, 1)) / omega,
        counts=counts,
        samples=samples,
    )


def curve_cross_section(curve: DeflectionCurve, edges: Array) -> Array:
    """偏向角の曲線から求めた、区間ごとの ``dσ/dΩ`` の理論値。

    区間に入る ρ の範囲の、円環の面積 ``π (ρ_{i+1}² − ρ_i²)`` を足し合わせる。
    曲線を細かく取り直して直線で補間するので、``χ(ρ)`` が単調でなくてもよい。
    """
    fine = np.linspace(curve.rho[0], curve.rho[-1], 400_000)
    middle = 0.5 * (fine[:-1] + fine[1:])
    chi = np.interp(middle, curve.rho, curve.chi)
    ring = np.pi * np.diff(fine**2)
    area, _ = np.histogram(chi, bins=edges, weights=ring)
    return area / solid_angles(edges)


def branch_cross_section(curve: DeflectionCurve, chi: Array) -> Array:
    """点ごとの ``dσ/dΩ = Σ_枝 (ρ / sin χ) / |dχ/dρ|``。"""
    values = curve.chi
    slope = np.gradient(values, curve.rho)
    result = np.zeros_like(chi, dtype=float)
    for j, target in enumerate(chi):
        diff = values - target
        for i in np.nonzero(diff[:-1] * diff[1:] < 0)[0]:
            w = diff[i] / (diff[i] - diff[i + 1])
            rho = curve.rho[i] + w * (curve.rho[i + 1] - curve.rho[i])
            derivative = slope[i] + w * (slope[i + 1] - slope[i])
            result[j] += rho / (np.sin(target) * abs(derivative))
    return result


def hard_sphere_angle_from_reflection(rho: Array, radius: float) -> Array:
    """球面で鏡のように跳ね返した速度の向きから散乱角を測る（公式は使わない）。"""
    vx, vy = hard_sphere_reflect(rho, radius)
    return np.abs(np.arctan2(vy, vx))


@dataclass(frozen=True)
class LabCrossSection:
    """実験室系で数えた、剛体球どうしの衝突の ``dσ/dΩ1``。"""

    mass_ratio: float  # m1 / m2
    counted: CrossSection
    expected: Array


def hard_sphere_lab_cross_section(
    m1: float,
    m2: float,
    radius: float = 1.0,
    n_bins: int = 30,
    samples: int = 200_000,
    seed: int = 1,
) -> LabCrossSection:
    """重心系で剛体球の散乱を起こし、実験室系の入射粒子の向きを数える。

    理論値は、崩壊の角度分布（``k = m1/m2``）に全断面積 ``π a²`` を掛けたもの。
    """
    k = m1 / m2
    limit = max_lab_angle(k)
    edges = np.linspace(0.0, limit, n_bins + 1)
    edges[-1] = min(limit + 1e-9, np.pi)

    def lab_angle_of(rho: Array) -> Array:
        vx, vy = hard_sphere_reflect(rho, radius)  # 重心系の向き
        direction = np.column_stack([vx, vy])
        u1, _ = final_velocities(m1, m2, [1.0, 0.0], [0.0, 0.0], direction)
        return np.arctan2(u1[:, 1], u1[:, 0])

    counted = count_cross_section(lab_angle_of, radius, edges, samples=samples, seed=seed)
    expected = np.pi * radius**2 * angle_bin_probability(edges, k) / solid_angles(edges)
    return LabCrossSection(k, counted, expected)
