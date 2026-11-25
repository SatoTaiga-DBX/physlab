"""慣性中心（重心）の実験: ばねでつないだ 3 質点を、回しながら一様重力の中へ投げる。

設計書: docs/design/006-conservation-laws.md ／ 教科書: 第 4 章 §慣性中心
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from physlab.application.mechanics.symmetry import run_cluster
from physlab.domain.mechanics.cluster import (
    ClusterParams,
    angular_momentum,
    as_vectors,
    center_of_mass,
    external_energy,
    kinetic_energy,
    momentum,
    spring_energy,
)
from physlab.numerics.integrators import Array


@dataclass(frozen=True)
class ThrowResult:
    params: ClusterParams
    t: Array
    q: Array
    v: Array
    center: Array  # 重心 R(t)  (len(t), 2)
    center_exact: Array  # R0 + V0 t + ½ g t²
    kinetic_center: Array  # ½ M V²（重心の運動）
    kinetic_relative: Array  # 重心に対する運動エネルギー
    spring: Array  # ばねのエネルギー
    gravity: Array  # 重力のポテンシャル M g Y
    internal_energy: Array  # 相対運動エネルギー + ばね（保存する）
    total_energy: Array
    angular_total: Array  # 原点まわり
    angular_center: Array  # R × P（重心の分）
    angular_relative: Array  # 重心まわり（保存する）


def initial_triangle(params: ClusterParams, stretch: float, spin: float, velocity: Array):
    """重心が原点の正三角形（辺 = 自然長 × ``stretch``）に、重心の速度と回転を与える。"""
    side = params.rest_length * stretch
    corners = side * np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3) / 2]])
    corners -= np.average(corners, axis=0, weights=params.masses)
    spin_velocity = spin * np.stack([-corners[:, 1], corners[:, 0]], axis=1)  # ω ẑ × r
    return corners.ravel(), (np.asarray(velocity, dtype=float) + spin_velocity).ravel()


def throw_cluster(
    velocity: tuple[float, float] = (3.0, 8.0),
    spin: float = 6.0,
    stretch: float = 1.15,
    stiffness: float = 200.0,
    g: float = 9.81,
    duration: float | None = None,
    dt: float = 1e-3,
) -> ThrowResult:
    """``duration`` を省略すると、重心が投げた高さに戻るまで（``2 V_y / g``）。"""
    params = ClusterParams(field="gravity", g=g, stiffness=stiffness)
    duration = 2.0 * velocity[1] / g if duration is None else duration
    q0, v0 = initial_triangle(params, stretch, spin, velocity)
    run = run_cluster(params, q0, v0, duration=duration, dt=dt)

    center = center_of_mass(run.q, params)
    total_momentum = momentum(run.v, params)
    center_velocity = total_momentum / params.total_mass
    relative = (as_vectors(run.v) - center_velocity[:, None, :]).reshape(run.v.shape)
    kinetic_center = 0.5 * params.total_mass * np.sum(center_velocity**2, axis=1)
    kinetic_relative = kinetic_energy(relative, params)
    spring = spring_energy(run.q, params)
    gravity = external_energy(run.t, run.q, params)
    angular_total = angular_momentum(run.q, run.v, params)
    angular_center = center[:, 0] * total_momentum[:, 1] - center[:, 1] * total_momentum[:, 0]
    relative_q = (as_vectors(run.q) - center[:, None, :]).reshape(run.q.shape)

    return ThrowResult(
        params=params,
        t=run.t,
        q=run.q,
        v=run.v,
        center=center,
        center_exact=center[0]
        + center_velocity[0] * run.t[:, None]
        + 0.5 * np.array([0.0, -g]) * run.t[:, None] ** 2,
        kinetic_center=kinetic_center,
        kinetic_relative=kinetic_relative,
        spring=spring,
        gravity=gravity,
        internal_energy=kinetic_relative + spring,
        total_energy=kinetic_center + kinetic_relative + spring + gravity,
        angular_total=angular_total,
        angular_center=angular_center,
        angular_relative=angular_momentum(relative_q, relative, params),
    )
