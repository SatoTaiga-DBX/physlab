"""ばねでつないだ質点の集まり（2 次元）と、外場による対称性の破れ。

設計書: docs/design/006-conservation-laws.md ／ 教科書: 第 4 章「保存法則」

内力は全ペアを結ぶばね ``U_int = Σ ½ k (|r_a − r_b| − ℓ)²``（相対位置だけに依る中心力）。
外場 ``field`` を切り替えると、壊れる対称性が変わる:

====================  ==========================================  ============================
field                 外場 U_ext                                   保存する量
====================  ==========================================  ============================
``none``              0                                            E, Px, Py, L
``gravity``           Σ m g y                                      E, Px
``trap``              ½ κ Σ m r²                                   E, L
``pulsing_trap``      ½ κ(t) Σ m r²,  κ = κ0 (1 + β sin Ω t)      L
====================  ==========================================  ============================

状態の並べ方: 位置 ``q = [x1, y1, x2, y2, …]``、速度 ``v`` も同じ並び（最後の軸が 2N 成分）。
時刻の配列に対してまとめて評価できる（先頭の軸が時刻）。
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]

FIELDS = ("none", "gravity", "trap", "pulsing_trap")


@dataclass(frozen=True)
class ClusterParams:
    masses: tuple[float, ...] = (1.0, 2.0, 1.5)
    stiffness: float = 20.0  # ばね定数 k
    rest_length: float = 1.0  # 自然長 ℓ
    field: str = "none"  # FIELDS のいずれか
    g: float = 1.0  # 一様重力（gravity）
    trap: float = 1.0  # 中心力の外場の強さ κ0（trap / pulsing_trap）
    pulse_depth: float = 0.5  # β
    pulse_frequency: float = 1.3  # Ω

    def __post_init__(self) -> None:
        if self.field not in FIELDS:
            raise ValueError(f"field は {FIELDS} のいずれか: {self.field!r}")

    @property
    def n(self) -> int:
        return len(self.masses)

    @property
    def total_mass(self) -> float:
        return float(sum(self.masses))


def as_vectors(x: Array) -> Array:
    """``[x1, y1, x2, y2, …]``（最後の軸が 2N）を ``(…, N, 2)`` にする。"""
    return x.reshape(*x.shape[:-1], -1, 2)


def _pairs(n: int) -> tuple[Array, Array]:
    return np.triu_indices(n, k=1)


def trap_strength(t: Array | float, p: ClusterParams) -> Array | float:
    """中心力の外場の強さ κ(t)。外場がそれ以外なら 0。"""
    if p.field == "trap":
        return p.trap + 0.0 * np.asarray(t)
    if p.field == "pulsing_trap":
        return p.trap * (1.0 + p.pulse_depth * np.sin(p.pulse_frequency * np.asarray(t)))
    return 0.0 * np.asarray(t)


def trap_strength_rate(t: Array | float, p: ClusterParams) -> Array | float:
    """dκ/dt。時間に依存するのは ``pulsing_trap`` だけ。"""
    if p.field == "pulsing_trap":
        return (
            p.trap * p.pulse_depth * p.pulse_frequency * np.cos(p.pulse_frequency * np.asarray(t))
        )
    return 0.0 * np.asarray(t)


def forces(t: float, q: Array, p: ClusterParams) -> tuple[Array, Array]:
    """``(内力, 外力)``。どちらも形状 ``(N, 2)``。"""
    r = as_vectors(q)
    m = np.asarray(p.masses)
    internal = np.zeros_like(r)
    i, j = _pairs(p.n)
    d = r[i] - r[j]
    length = np.linalg.norm(d, axis=-1, keepdims=True)
    pull = -p.stiffness * (length - p.rest_length) * d / length  # i が j から受ける力
    np.add.at(internal, i, pull)
    np.add.at(internal, j, -pull)

    external = np.zeros_like(r)
    if p.field == "gravity":
        external[:, 1] = -m * p.g
    elif p.field in ("trap", "pulsing_trap"):
        external = -float(trap_strength(t, p)) * m[:, None] * r
    return internal, external


def accel(p: ClusterParams) -> Callable[[float, Array], Array]:
    """``a(t, q)``。速度ベルレ法に渡す。"""
    m = np.repeat(np.asarray(p.masses), 2)

    def a(t: float, q: Array) -> Array:
        internal, external = forces(t, q, p)
        return (internal + external).ravel() / m

    return a


def kinetic_energy(v: Array, p: ClusterParams) -> Array:
    return 0.5 * np.sum(np.asarray(p.masses) * np.sum(as_vectors(v) ** 2, axis=-1), axis=-1)


def spring_energy(q: Array, p: ClusterParams) -> Array:
    """内部のポテンシャル U_int。"""
    r = as_vectors(q)
    i, j = _pairs(p.n)
    length = np.linalg.norm(r[..., i, :] - r[..., j, :], axis=-1)
    return 0.5 * p.stiffness * np.sum((length - p.rest_length) ** 2, axis=-1)


def external_energy(t: Array | float, q: Array, p: ClusterParams) -> Array:
    """外場のポテンシャル U_ext。"""
    r = as_vectors(q)
    m = np.asarray(p.masses)
    if p.field == "gravity":
        return p.g * np.sum(m * r[..., 1], axis=-1)
    squared = np.sum(m * np.sum(r**2, axis=-1), axis=-1)
    return 0.5 * trap_strength(t, p) * squared


def energy(t: Array | float, q: Array, v: Array, p: ClusterParams) -> Array:
    """``E = T + U_int + U_ext``。"""
    return kinetic_energy(v, p) + spring_energy(q, p) + external_energy(t, q, p)


def momentum(v: Array, p: ClusterParams) -> Array:
    """全運動量 ``P = Σ m v``（形状 ``(…, 2)``）。"""
    return np.sum(np.asarray(p.masses)[:, None] * as_vectors(v), axis=-2)


def angular_momentum(q: Array, v: Array, p: ClusterParams) -> Array:
    """原点まわりの角運動量の z 成分 ``L = Σ m (x v_y − y v_x)``。"""
    r, u = as_vectors(q), as_vectors(v)
    cross = r[..., 0] * u[..., 1] - r[..., 1] * u[..., 0]
    return np.sum(np.asarray(p.masses) * cross, axis=-1)


def center_of_mass(q: Array, p: ClusterParams) -> Array:
    """重心 ``R = Σ m r / M``（形状 ``(…, 2)``）。"""
    return np.sum(np.asarray(p.masses)[:, None] * as_vectors(q), axis=-2) / p.total_mass


def predicted_rates(t: Array, q: Array, p: ClusterParams) -> tuple[Array, Array, Array]:
    """対称性の破れから決まる ``(dE/dt, dP/dt, dL/dt)``。内力は寄与しない。

    * ``dE/dt = ∂U_ext/∂t``（時間の一様性の破れ）
    * ``dP/dt = Σ F_ext``（並進対称性の破れ）
    * ``dL/dt = Σ (r × F_ext)_z``（回転対称性の破れ）
    """
    t = np.asarray(t, dtype=float)
    r = as_vectors(q)
    m = np.asarray(p.masses)
    if p.field == "gravity":
        external = np.zeros_like(r)
        external[..., 1] = -m * p.g
    else:
        external = -np.asarray(trap_strength(t, p))[..., None, None] * m[:, None] * r
    d_energy = 0.5 * trap_strength_rate(t, p) * np.sum(m * np.sum(r**2, axis=-1), axis=-1)
    d_momentum = np.sum(external, axis=-2)
    d_angular = np.sum(r[..., 0] * external[..., 1] - r[..., 1] * external[..., 0], axis=-1)
    return d_energy, d_momentum, d_angular
