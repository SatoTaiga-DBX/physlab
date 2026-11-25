"""面積速度の実験: 中心力の軌道を、等しい時間の扇形に区切って面積と弧の長さを比べる。

設計書: docs/design/006-conservation-laws.md ／ 教科書: 第 4 章 §角運動量
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from physlab.domain.mechanics.central import (
    CentralParams,
    accel,
    angular_momentum,
    energy,
)
from physlab.numerics.integrators import Array, velocity_verlet


@dataclass(frozen=True)
class CentralOrbit:
    params: CentralParams
    t: Array
    q: Array  # (len(t), 2)
    v: Array
    angular: Array
    energy: Array
    step_area: Array  # 1 ステップで動径が掃く三角形の面積（len(t) − 1）


@dataclass(frozen=True)
class Sectors:
    """軌道を等しい時間の区間に分けた扇形。"""

    orbit: CentralOrbit
    edges: Array  # 区切りの添字（len = n_sectors + 1）
    areas: Array
    arc_lengths: Array
    duration: float  # 1 区間の時間


def central_orbit(
    params: CentralParams | None = None,
    q0: tuple[float, float] = (1.0, 0.0),
    v0: tuple[float, float] = (0.0, 0.9),
    duration: float = 6.0,
    dt: float = 5e-4,
) -> CentralOrbit:
    params = CentralParams() if params is None else params
    t = np.arange(0.0, duration + dt / 2, dt)
    q, v = velocity_verlet(accel(params), q0, v0, t)
    cross = q[:-1, 0] * q[1:, 1] - q[:-1, 1] * q[1:, 0]
    return CentralOrbit(
        params=params,
        t=t,
        q=q,
        v=v,
        angular=angular_momentum(q, v, params),
        energy=energy(q, v, params),
        step_area=0.5 * np.abs(cross),
    )


def time_of_one_revolution(orbit: CentralOrbit) -> float:
    """中心のまわりを初めてちょうど 1 周（偏角が 2π 進む）するまでの時間。"""
    angle = np.unwrap(np.arctan2(orbit.q[:, 1], orbit.q[:, 0]))
    turned = np.abs(angle - angle[0])
    if turned[-1] < 2 * np.pi:
        raise ValueError("軌道が 1 周していない。duration を長くする")
    return float(orbit.t[np.searchsorted(turned, 2 * np.pi)])


def equal_time_sectors(
    orbit: CentralOrbit, n_sectors: int = 12, span: float | None = None
) -> Sectors:
    """最初の ``span``（省略で全体）を ``n_sectors`` 等分し、各区間の面積と弧の長さを求める。"""
    dt = float(orbit.t[1] - orbit.t[0])
    total = orbit.t.size - 1 if span is None else int(round(span / dt))
    steps = total // n_sectors
    edges = np.arange(n_sectors + 1) * steps
    segment = np.linalg.norm(np.diff(orbit.q, axis=0), axis=1)
    areas = np.array(
        [orbit.step_area[a:b].sum() for a, b in zip(edges[:-1], edges[1:], strict=True)]
    )
    arcs = np.array([segment[a:b].sum() for a, b in zip(edges[:-1], edges[1:], strict=True)])
    return Sectors(orbit, edges, areas, arcs, steps * dt)
