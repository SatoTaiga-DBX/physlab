"""中心力による偏向角を 2 通りで求める: 軌道の式の求積と、運動方程式の時間積分。

設計書: docs/design/007-collisions.md ／ 教科書: 第 5 章 §散乱

* 求積: ``φ0 = ∫0^{u_max} ρ du / √F(u)``、``Θ = π − 2 φ0``。
* 時間積分: 有限の距離 ``R0`` から出発する。無限遠の衝突パラメータが正しく ρ になるように、
  ``R0`` までに回る角を同じ求積で求めて初期条件を置く。終わりも無限遠まで角を延ばす。

入射は左（x = −∞）から +x 向きで、``y = +ρ`` の側を通る。偏向角 Θ は反時計回りが正
（斥力で正、引力で負）。

仮定: ``F(u)`` の正の零点は 1 つ（単調な斥力、またはクーロン型の引力）。
引力でほぼ正面から入る軌道（クーロン引力で ``ρ`` が最接近の目安 ``d`` の 1/100 以下）は、
中心の近くで鋭く折り返すので、求積の精度が落ちる。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from physlab.domain.mechanics.potentials import CentralPotential
from physlab.domain.mechanics.scattering import accel, radial_function
from physlab.numerics.integrators import Array, velocity_verlet
from physlab.numerics.quadrature import integrate_inverse_sqrt_end, integrate_smooth
from physlab.numerics.roots import bisect


@dataclass(frozen=True)
class Orbits:
    """固定した中心のまわりを通る粒子の束。"""

    potential: CentralPotential
    mass: float
    speed: float
    rho: Array  # (M,)
    t: Array
    x: Array  # (len(t), M)
    y: Array
    vx: Array
    vy: Array


def turning_u(rho: ArrayLike, potential: CentralPotential, energy: float) -> Array:
    """最接近の ``u_max = 1/r_min``（``F(u) = 0`` の零点）。``ρ > 0`` の配列。"""
    rho = np.atleast_1d(np.asarray(rho, dtype=float))
    if np.any(rho <= 0):
        raise ValueError("rho は正の値（正面衝突は ρ → 0 の極限）")

    def f(u: Array) -> Array:
        return radial_function(u, rho, potential, energy)

    upper = 1.0 / rho
    for _ in range(200):  # 引力では 1/ρ の先に零点があるので、挟めるまで広げる
        outside = f(upper) > 0
        if not outside.any():
            break
        upper = np.where(outside, 2.0 * upper, upper)
    return bisect(f, np.zeros_like(rho), upper)


def deflection(rho: ArrayLike, potential: CentralPotential, energy: float, n: int = 120) -> Array:
    """偏向角 Θ（符号つき、斥力で正）。"""
    rho = np.atleast_1d(np.asarray(rho, dtype=float))
    u_max = turning_u(rho, potential, energy)
    column = rho[:, None]

    def integrand(u: Array) -> Array:
        return column / np.sqrt(radial_function(u, column, potential, energy))

    # 式の切り替わる半径があればそこで、なければ中ほどで区間を分ける
    split = 0.5 * u_max
    if potential.breakpoint_radius is not None:
        edge = 1.0 / potential.breakpoint_radius
        split = np.where(edge < u_max, edge, split)
    # 前半も端に節点を集める。球の表面をかすめる軌道では、分け目で F がほとんど 0 になるため
    phi = integrate_inverse_sqrt_end(integrand, np.zeros_like(u_max), split, n=n)
    phi += integrate_inverse_sqrt_end(integrand, split, u_max, n=n)
    return np.pi - 2.0 * phi


def deflection_in_chunks(
    rho: ArrayLike,
    potential: CentralPotential,
    energy: float,
    n: int = 64,
    chunk: int = 20_000,
) -> Array:
    """多数の ρ の偏向角。メモリを使いすぎないように分けて計算する。"""
    rho = np.atleast_1d(np.asarray(rho, dtype=float))
    parts = [
        deflection(rho[i : i + chunk], potential, energy, n=n) for i in range(0, rho.size, chunk)
    ]
    return np.concatenate(parts)


def swept_angle(
    rho: ArrayLike, radius: ArrayLike, potential: CentralPotential, energy: float
) -> Array:
    """半径 ``radius`` から無限遠までに、中心のまわりを回る角 ``∫0^{1/R} ρ du/√F``。"""
    rho = np.atleast_1d(np.asarray(rho, dtype=float))
    radius = np.broadcast_to(np.asarray(radius, dtype=float), rho.shape)
    column = rho[:, None]

    def integrand(u: Array) -> Array:
        return column / np.sqrt(radial_function(u, column, potential, energy))

    return integrate_smooth(integrand, np.zeros_like(rho), 1.0 / radius, n=80)


def incoming_state(
    potential: CentralPotential,
    mass: float,
    speed: float,
    rho: ArrayLike,
    start_radius: ArrayLike,
) -> tuple[Array, Array]:
    """距離 ``start_radius`` の入射側の初期条件 ``(q0, v0)``（``[x…, y…]`` の並び）。

    ``start_radius`` は粒子ごとの配列でもよい。
    """
    rho = np.atleast_1d(np.asarray(rho, dtype=float))
    energy = 0.5 * mass * speed**2
    angle = np.pi - swept_angle(rho, start_radius, potential, energy)
    cos_a, sin_a = np.cos(angle), np.sin(angle)
    u = np.full_like(rho, 1.0 / start_radius)
    radial = -speed * np.sqrt(radial_function(u, rho, potential, energy))
    around = -speed * rho / start_radius  # 時計回り（角運動量 −m v ρ）
    q0 = np.concatenate([start_radius * cos_a, start_radius * sin_a])
    v0 = np.concatenate([radial * cos_a - around * sin_a, radial * sin_a + around * cos_a])
    return q0, v0


def trace_orbits(
    potential: CentralPotential,
    mass: float,
    speed: float,
    rho: ArrayLike,
    start_radius: ArrayLike = 12.0,
    duration: float | None = None,
    dt: float = 1e-3,
) -> Orbits:
    """粒子の束の軌道を速度ベルレ法で解く。``ρ`` が負なら中心の下側を通る。

    ``start_radius`` を ``√(X² + ρ²)`` の配列にすると、束の先頭がほぼ ``x = −X`` にそろう。
    """
    rho = np.atleast_1d(np.asarray(rho, dtype=float))
    side = np.where(rho < 0, -1.0, 1.0)  # 下側は上側の鏡像
    q0, v0 = incoming_state(potential, mass, speed, np.abs(rho), start_radius)
    q0[rho.size :] *= side
    v0[rho.size :] *= side
    duration = 2.2 * float(np.max(start_radius)) / speed if duration is None else duration
    t = np.arange(0.0, duration + dt / 2, dt)
    q, v = velocity_verlet(accel(mass, potential), q0, v0, t)
    n = rho.size
    return Orbits(potential, mass, speed, rho, t, q[:, :n], q[:, n:], v[:, :n], v[:, n:])


def orbit_deflection(orbits: Orbits) -> Array:
    """時間積分の結果から偏向角を測る（最後の点から無限遠までの角を足す）。

    Raises:
        ValueError: まだ遠ざかり始めていない粒子がある（積分時間が足りない）。
    """
    x, y = orbits.x, orbits.y
    radial_speed = x[-1] * orbits.vx[-1] + y[-1] * orbits.vy[-1]
    if np.any(radial_speed <= 0):
        raise ValueError("遠ざかっていない粒子がある。duration を長くする")
    side = np.where(orbits.rho < 0, -1.0, 1.0)
    angle = np.unwrap(np.arctan2(side * y, x), axis=0)[-1]
    energy = 0.5 * orbits.mass * orbits.speed**2
    tail = swept_angle(np.abs(orbits.rho), np.hypot(x[-1], y[-1]), orbits.potential, energy)
    return side * (angle - tail)


def closest_distance(orbits: Orbits) -> Array:
    """軌道が中心に最も近づいた距離（記録した点の最小値）。"""
    return np.hypot(orbits.x, orbits.y).min(axis=0)
