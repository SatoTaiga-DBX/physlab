"""ガリレイの相対性のユースケース: 電車の中の投げ上げ、作用のずれ、二体問題を 2 つの慣性系で。

設計書: docs/design/005-analytical-mechanics-ch1.md ／ 教科書: 第 3 章 §ガリレイの相対性原理
返すのは dataclass（配列と派生量）だけ。これを描画するのは infrastructure の仕事。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from physlab.domain.mechanics.galilean import free_particle_action_shift
from physlab.domain.mechanics.lagrangian import action, free_particle, sine_mode
from physlab.domain.mechanics.two_body import (
    TwoBodyParams,
    accel,
    center_of_mass,
    energy,
    initial_state_cm,
    total_momentum,
)
from physlab.numerics.integrators import Array, integrate, velocity_verlet

# ---- 電車の中の投げ上げ ------------------------------------------------------------------


@dataclass(frozen=True)
class TrainThrow:
    """電車の中で球を真上に投げ上げる。電車系 K' と地面系 K で同じ運動を記述する。"""

    train_speed: float  # K' の K に対する速度（x 方向）[m/s]
    throw_speed: float  # K' での鉛直上向きの初速 [m/s]
    g: float
    t: Array
    train: Array  # K' での位置 (len(t), 2)
    ground: Array  # K での位置 (len(t), 2)
    accel_train: Array  # 位置の 2 階差分から求めた加速度（内部の点）
    accel_ground: Array
    kinetic_train: Array  # 単位質量あたりの運動エネルギー
    kinetic_ground: Array


def throw_in_train(
    train_speed: float,
    throw_speed: float,
    g: float = 9.81,
    duration: float | None = None,
    dt: float = 1e-3,
) -> TrainThrow:
    """ニュートンの運動方程式（``ẍ = 0``, ``ÿ = −g``）を両方の系で独立に RK4 で解く。"""
    duration = 2.0 * throw_speed / g if duration is None else duration
    t = np.arange(0.0, duration + dt / 2, dt)

    def f(_t: float, y: Array) -> Array:  # y = [x, y, vx, vy]。どちらの系でも同じ式
        return np.array([y[2], y[3], 0.0, -g])

    train = integrate(f, [0.0, 0.0, 0.0, throw_speed], t)  # K': まっすぐ上へ
    ground = integrate(f, [0.0, 0.0, train_speed, throw_speed], t)  # K: 横にも動いている
    return TrainThrow(
        train_speed=train_speed,
        throw_speed=throw_speed,
        g=g,
        t=t,
        train=train[:, :2],
        ground=ground[:, :2],
        accel_train=np.diff(train[:, :2], 2, axis=0) / dt**2,
        accel_ground=np.diff(ground[:, :2], 2, axis=0) / dt**2,
        kinetic_train=0.5 * np.sum(train[:, 2:] ** 2, axis=1),
        kinetic_ground=0.5 * np.sum(ground[:, 2:] ** 2, axis=1),
    )


# ---- 作用のずれ（自由粒子） ---------------------------------------------------------------


@dataclass(frozen=True)
class ActionShift:
    """K での作用 ``S`` と、K' での作用 ``S'``。差は経路によらない定数。"""

    mass: float
    frame_speed: float
    t: Array
    eps: Array
    action_k: Array  # S(ε)
    action_kprime: Array  # S'(ε)
    shift: Array  # S − S'
    shift_theory: float  # m V Δq' + ½ m V² τ


def action_shift_under_boost(
    mass: float,
    frame_speed: float,
    t1: float,
    q1: float,
    t2: float,
    q2: float,
    eps: Array,
    n_points: int = 4001,
) -> ActionShift:
    """両端を固定した経路の族 ``q = (直線) + ε sin(π s/τ)`` で、``S`` と ``S'`` を比べる。

    ``S'`` は K' での経路 ``q' = q − V t``、``q̇' = q̇ − V`` に対する作用。
    """
    system = free_particle(mass)
    t = np.linspace(t1, t2, n_points)
    tau = t2 - t1
    eta = sine_mode(t, t1, t2, 1)
    eta_dot = np.pi / tau * np.cos(np.pi * (t - t1) / tau)
    line = q1 + (q2 - q1) * (t - t1) / tau
    line_dot = (q2 - q1) / tau
    eps = np.asarray(eps, dtype=float)

    s_k = np.empty_like(eps)
    s_kp = np.empty_like(eps)
    for i, e in enumerate(eps):
        q, qdot = line + e * eta, line_dot + e * eta_dot
        s_k[i] = action(system, t, q, qdot)
        s_kp[i] = action(system, t, q - frame_speed * t, qdot - frame_speed)
    dq_prime = (q2 - frame_speed * t2) - (q1 - frame_speed * t1)
    return ActionShift(
        mass=mass,
        frame_speed=frame_speed,
        t=t,
        eps=eps,
        action_k=s_k,
        action_kprime=s_kp,
        shift=s_k - s_kp,
        shift_theory=free_particle_action_shift(mass, frame_speed, dq_prime, tau),
    )


# ---- 二体問題を 2 つの慣性系で -------------------------------------------------------------


@dataclass(frozen=True)
class TwoBodyFrames:
    """重心系 K' と、それに対して速度 V で動く実験室系 K で、同じ二体運動を独立に解いた結果。"""

    params: TwoBodyParams
    frame_velocity: Array  # K' の K に対する速度 (2,)
    t: Array
    q_cm: Array  # K' での位置 (len(t), 4)
    v_cm: Array
    q_lab: Array  # K での位置（K で独立に積分）
    v_lab: Array
    frame_difference: Array  # max |q_lab − (q_cm + V t)| の時間変化
    energy_cm: Array
    energy_lab: Array
    momentum_lab: Array  # K での全運動量 (len(t), 2)
    center_lab: Array  # K での重心 (len(t), 2)


def simulate_two_body_frames(
    params: TwoBodyParams, frame_velocity: Array, duration: float, dt: float = 1e-3
) -> TwoBodyFrames:
    """重心系 K' で解き、実験室系 K では「初速に V を足しただけ」の初期条件で別に解く。

    速度ベルレ法は加速度が相対位置だけに依存する系で V を足しても更新式が形を変えないので、
    2 つの解は丸め誤差の範囲で ``q_lab = q_cm + V t`` に一致する。
    """
    v_frame = np.asarray(frame_velocity, dtype=float)
    t = np.arange(0.0, duration + dt / 2, dt)
    acc = accel(params)
    q0, v0 = initial_state_cm(params)
    q_cm, v_cm = velocity_verlet(acc, q0, v0, t)
    q_lab, v_lab = velocity_verlet(acc, q0, v0 + np.tile(v_frame, 2), t)

    expected = q_cm + t[:, None] * np.tile(v_frame, 2)[None, :]
    return TwoBodyFrames(
        params=params,
        frame_velocity=v_frame,
        t=t,
        q_cm=q_cm,
        v_cm=v_cm,
        q_lab=q_lab,
        v_lab=v_lab,
        frame_difference=np.max(np.abs(q_lab - expected), axis=1),
        energy_cm=energy(q_cm, v_cm, params),
        energy_lab=energy(q_lab, v_lab, params),
        momentum_lab=total_momentum(v_lab, params),
        center_lab=center_of_mass(q_lab, params),
    )
