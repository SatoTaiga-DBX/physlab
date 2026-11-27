"""固定した力の中心による散乱（相対運動）の物理と、解析解のある場合の公式。

設計書: docs/design/007-collisions.md ／ 教科書: 第 5 章 §散乱・§ラザフォードの公式・§微小角

記号: 換算質量 ``mass``、無限遠の速さ ``speed``、エネルギー ``E = ½ mass speed²``、
衝突パラメータ ``ρ``、偏向角 ``Θ``（斥力で正）、散乱角 ``χ = |Θ|``。

軌道の形は ``u = 1/r`` で表すと ``dφ = ρ du / √F(u)``、
``F(u) = 1 − ρ² u² − U(1/u)/E``。最接近は ``F = 0`` で、そこまでの角 φ0 を使って
``Θ = π − 2 φ0``。

状態の並べ方（粒子の束を同時に解くとき）: 位置 ``q = [x1…xN, y1…yN]``。
"""

from __future__ import annotations

from collections.abc import Callable
from math import gamma, pi, sqrt

import numpy as np
from numpy.typing import ArrayLike, NDArray

from physlab.domain.mechanics.potentials import CentralPotential

Array = NDArray[np.float64]


def kinetic_at_infinity(mass: float, speed: float) -> float:
    return 0.5 * mass * speed**2


def radial_function(u: Array, rho: Array, potential: CentralPotential, energy: float) -> Array:
    """``F(u) = 1 − ρ² u² − U(1/u)/E``。``(dr/dt)² = v∞² F`` なので、``F ≥ 0`` の範囲を動く。"""
    return 1.0 - (rho * u) ** 2 - potential.value(1.0 / u) / energy


def accel(mass: float, potential: CentralPotential) -> Callable[[float, Array], Array]:
    """固定した中心のまわりの粒子の束 ``q = [x…, y…]`` の加速度。"""

    def a(_t: float, q: Array) -> Array:
        n = q.shape[-1] // 2
        x, y = q[..., :n], q[..., n:]
        r = np.hypot(x, y)
        radial = -potential.derivative(r) / (mass * r)
        return np.concatenate([radial * x, radial * y], axis=-1)

    return a


# ---------------------------------------------------------------- クーロン力（ラザフォード）


def head_on_distance(alpha: float, mass: float, speed: float) -> float:
    """``d = |α| / E``。斥力で正面から来た粒子が止まって引き返す距離。"""
    return abs(alpha) / kinetic_at_infinity(mass, speed)


def rutherford_angle(rho: ArrayLike, d: float) -> Array:
    """散乱角 ``χ = 2 arctan(d / 2ρ)``。引力と斥力で同じ大きさ。"""
    return 2.0 * np.arctan2(d, 2.0 * np.asarray(rho, dtype=float))


def rutherford_impact_parameter(chi: ArrayLike, d: float) -> Array:
    """``ρ = (d/2) cot(χ/2)``。"""
    return 0.5 * d / np.tan(0.5 * np.asarray(chi, dtype=float))


def rutherford_cross_section(chi: ArrayLike, d: float) -> Array:
    """``dσ/dΩ = (d/4)² / sin⁴(χ/2)``。"""
    return (0.25 * d) ** 2 / np.sin(0.5 * np.asarray(chi, dtype=float)) ** 4


def coulomb_closest_approach(chi: ArrayLike, d: float) -> Array:
    """斥力の最接近距離 ``r_min = (d/2) (1 + 1/sin(χ/2))``。"""
    return 0.5 * d * (1.0 + 1.0 / np.sin(0.5 * np.asarray(chi, dtype=float)))


def max_energy_transfer(m1: float, m2: float, speed: float) -> float:
    """静止した標的に渡せる最大のエネルギー ``2 μ² v² / m2``。"""
    mu = m1 * m2 / (m1 + m2)
    return 2.0 * mu**2 * speed**2 / m2


def energy_transfer_cross_section(
    epsilon: ArrayLike, alpha: float, m2: float, speed: float
) -> Array:
    """クーロン散乱で標的に ``ε`` を渡す断面積 ``dσ/dε = 2π α² / (m2 v² ε²)``。"""
    return 2.0 * pi * alpha**2 / (m2 * speed**2 * np.asarray(epsilon, dtype=float) ** 2)


def ball_graze_angle(radius_over_d: float) -> float:
    """一様に帯電した球（``R > d``）で、点電荷と同じに見える散乱角の上限 ``χ_g``。

    最接近がちょうど表面（``r_min = R``）になる角で、``sin(χ_g/2) = 1/(2R/d − 1)``。
    ``R ≤ d`` ならどの軌道も球に入らないので π。
    """
    if radius_over_d <= 1.0:
        return float(pi)
    return float(2.0 * np.arcsin(1.0 / (2.0 * radius_over_d - 1.0)))


# ---------------------------------------------------------------- 逆 2 乗ポテンシャル α/r²


def inverse_square_deflection(rho: ArrayLike, alpha: float, energy: float) -> Array:
    """``U = α/r²``（``α > 0``）の偏向角 ``π (1 − ρ / √(ρ² + α/E))``。厳密解の照合用。"""
    rho = np.asarray(rho, dtype=float)
    return pi * (1.0 - rho / np.sqrt(rho**2 + alpha / energy))


# ---------------------------------------------------------------- 剛体球


def hard_sphere_angle(rho: ArrayLike, radius: float) -> Array:
    """半径 a の固い球で鏡のように跳ね返るときの散乱角 ``π − 2 arcsin(ρ/a)``（``ρ ≥ a`` は 0）。"""
    ratio = np.clip(np.asarray(rho, dtype=float) / radius, 0.0, 1.0)
    return pi - 2.0 * np.arcsin(ratio)


def hard_sphere_reflect(rho: ArrayLike, radius: float) -> tuple[Array, Array]:
    """向き (1, 0) で来た粒子が球面で跳ね返った後の向き ``v − 2 (v·n) n``。

    当たらない（``ρ ≥ a``）粒子はそのまま直進する。
    """
    rho = np.asarray(rho, dtype=float)
    nx = -np.sqrt(np.maximum(1.0 - (rho / radius) ** 2, 0.0))  # 当たった点の法線
    ny = rho / radius
    hit = np.abs(rho) < radius
    along = np.where(hit, nx, 0.0)  # v·n（v = (1, 0)）
    return 1.0 - 2.0 * along * nx, -2.0 * along * ny


# ---------------------------------------------------------------- 微小角（力積の近似）


def power_law_small_angle(rho: ArrayLike, alpha: float, n: float, mass: float, speed: float):
    """``U = α/rⁿ`` の微小角の偏向 ``θ = 2√π α Γ((n+1)/2) / (m v² ρⁿ Γ(n/2))``。"""
    factor = 2.0 * sqrt(pi) * gamma(0.5 * (n + 1.0)) / gamma(0.5 * n)
    return factor * alpha / (mass * speed**2 * np.asarray(rho, dtype=float) ** n)


def power_law_small_angle_cross_section(theta: ArrayLike, coefficient: float, n: float) -> Array:
    """``θ = A ρ^{−n}`` のときの微小角の ``dσ/dΩ = (1/n) A^{2/n} θ^{−2−2/n}``。"""
    theta = np.asarray(theta, dtype=float)
    return coefficient ** (2.0 / n) / n * theta ** (-2.0 - 2.0 / n)


def ball_small_angle(rho: ArrayLike, d: float, radius: float) -> Array:
    """一様に帯電した球の微小角の偏向 ``(d/ρ) f(ρ)``。

    ``f`` は、半径 ρ の円柱（入射方向に平行）の内側にある電荷の割合
    ``1 − (1 − ρ²/R²)^{3/2}``。球の外を通る（``ρ ≥ R``）と ``f = 1``。
    """
    rho = np.asarray(rho, dtype=float)
    inside = 1.0 - np.clip(1.0 - (rho / radius) ** 2, 0.0, None) ** 1.5
    return d * inside / rho
