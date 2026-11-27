"""弾性衝突の実験: 短距離の力で押し合う 2 粒子の運動方程式を解き、衝突後の運動量を調べる。

設計書: docs/design/007-collisions.md ／ 教科書: 第 5 章 §弾性衝突

標的（粒子 2）は原点に静止し、入射粒子（粒子 1）は左から、衝突パラメータ ρ だけ上にずれて来る。
同じ質量・同じ力で ρ の違う衝突は、まとめて 1 回の積分で解く。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from physlab.domain.mechanics.collision import accel, energy
from physlab.domain.mechanics.potentials import (
    CentralPotential,
    Gaussian,
    PowerLaw,
    SoftDisk,
    Yukawa,
)
from physlab.numerics.integrators import Array, velocity_verlet

# 形のまったく違う短距離の斥力（衝突後の運動量が同じ円に乗ることを見る）
DEMO_POTENTIALS: dict[str, CentralPotential] = {
    "soft-disk": SoftDisk(stiffness=4.0, radius=1.0),
    "gaussian": Gaussian(height=1.0, width=0.8),
    "yukawa": Yukawa(alpha=1.0, length=0.6),
    "steep-core": PowerLaw(alpha=0.1, n=12.0),
}


@dataclass(frozen=True)
class CollisionBatch:
    m1: float
    m2: float
    speed: float
    potential: CentralPotential
    impact: Array  # (M,)
    t: Array
    q: Array  # (len(t), M, 4)  [x1, y1, x2, y2]
    v: Array
    energy: Array  # (len(t), M)


@dataclass(frozen=True)
class CollisionOutcome:
    """衝突の前後だけから読み取った量。"""

    impact: Array
    chi: Array  # 重心系で相対速度が回った角（0〜π）
    theta1: Array  # 実験室系で入射粒子が曲がった角
    theta2: Array  # 標的が飛び出した角（力を受けなかったときは NaN）
    p1: Array  # (M, 2)
    p2: Array
    transfer: Array  # E2' / E1


def simulate_collisions(
    m1: float,
    m2: float,
    potential: CentralPotential,
    impact: ArrayLike,
    speed: float = 1.0,
    start: float = 8.0,
    duration: float = 20.0,
    dt: float = 2e-3,
) -> CollisionBatch:
    rho = np.atleast_1d(np.asarray(impact, dtype=float))
    count = rho.size
    q0 = np.zeros((count, 4))
    q0[:, 0] = -start
    q0[:, 1] = rho
    v0 = np.zeros((count, 4))
    v0[:, 0] = speed
    t = np.arange(0.0, duration + dt / 2, dt)
    q, v = velocity_verlet(accel(m1, m2, potential), q0.ravel(), v0.ravel(), t)
    return CollisionBatch(
        m1=m1,
        m2=m2,
        speed=speed,
        potential=potential,
        impact=rho,
        t=t,
        q=q.reshape(t.size, count, 4),
        v=v.reshape(t.size, count, 4),
        energy=energy(q, v, m1, m2, potential),
    )


def outcome(batch: CollisionBatch) -> CollisionOutcome:
    final = batch.v[-1]
    u1, u2 = final[:, 0:2], final[:, 2:4]
    relative = u1 - u2
    chi = np.abs(np.arctan2(relative[:, 1], relative[:, 0]))
    theta1 = np.abs(np.arctan2(u1[:, 1], u1[:, 0]))
    moved = np.linalg.norm(u2, axis=1) > 1e-9
    theta2 = np.where(moved, np.abs(np.arctan2(u2[:, 1], u2[:, 0])), np.nan)
    transfer = batch.m2 * np.sum(u2**2, axis=1) / (batch.m1 * batch.speed**2)
    return CollisionOutcome(
        impact=batch.impact,
        chi=chi,
        theta1=theta1,
        theta2=theta2,
        p1=batch.m1 * u1,
        p2=batch.m2 * u2,
        transfer=transfer,
    )


def circle_experiment(
    m1: float,
    m2: float,
    impact: ArrayLike,
    potentials: dict[str, CentralPotential] | None = None,
    **kwargs,
) -> dict[str, CollisionOutcome]:
    """力の形を変えて同じ衝突を行う。"""
    potentials = DEMO_POTENTIALS if potentials is None else potentials
    return {
        name: outcome(simulate_collisions(m1, m2, potential, impact, **kwargs))
        for name, potential in potentials.items()
    }
