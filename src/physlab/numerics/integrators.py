"""固定刻みの常微分方程式（ODE）ソルバ。

学習目的で「手法ごとの精度・保存則の守られ方の違い」を比べられるよう、
あえて自前で小さく実装している。本格的な計算には scipy.integrate.solve_ivp を使うこと。

2 つの流儀を用意している:

* 一般の 1 階 ODE  ``dy/dt = f(t, y)``  → :func:`integrate`（euler / rk4）
* 2 階の力学系     ``d²q/dt² = a(t, q)`` → :func:`velocity_verlet`（シンプレクティック）
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray

Array = NDArray[np.float64]
Deriv = Callable[[float, Array], Array]  # f(t, y) -> dy/dt
Accel = Callable[[float, Array], Array]  # a(t, q) -> d²q/dt²
Stepper = Callable[[Deriv, float, Array, float], Array]


def euler_step(f: Deriv, t: float, y: Array, dt: float) -> Array:
    """前進オイラー法（1 次精度）。"""
    return y + dt * f(t, y)


def rk4_step(f: Deriv, t: float, y: Array, dt: float) -> Array:
    """古典的 4 次ルンゲ＝クッタ法（4 次精度）。"""
    k1 = f(t, y)
    k2 = f(t + dt / 2, y + dt / 2 * k1)
    k3 = f(t + dt / 2, y + dt / 2 * k2)
    k4 = f(t + dt, y + dt * k3)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


STEPPERS: dict[str, Stepper] = {"euler": euler_step, "rk4": rk4_step}


def integrate(f: Deriv, y0: ArrayLike, t: ArrayLike, method: str = "rk4") -> Array:
    """``dy/dt = f(t, y)`` を時刻列 ``t`` 上で解く。

    Returns:
        形状 ``(len(t), len(y0))`` の配列。``t`` は等間隔でなくてもよい。
    """
    step = STEPPERS[method]
    t = np.asarray(t, dtype=float)
    y0 = np.atleast_1d(np.asarray(y0, dtype=float))
    y = np.empty((t.size, y0.size))
    y[0] = y0
    for i in range(t.size - 1):
        y[i + 1] = step(f, t[i], y[i], t[i + 1] - t[i])
    return y


def velocity_verlet(acc: Accel, q0: ArrayLike, v0: ArrayLike, t: ArrayLike) -> tuple[Array, Array]:
    """速度ベルレ法（2 次精度・シンプレクティック）。

    エネルギーが長時間にわたって有界に振動し、ドリフトしないのが特徴。
    ``acc`` が速度に依存しない（保存力のみ）場合に使える。

    Returns:
        ``(q, v)``。それぞれ形状 ``(len(t), len(q0))``。
    """
    t = np.asarray(t, dtype=float)
    q0 = np.atleast_1d(np.asarray(q0, dtype=float))
    v0 = np.atleast_1d(np.asarray(v0, dtype=float))
    q = np.empty((t.size, q0.size))
    v = np.empty_like(q)
    q[0], v[0] = q0, v0
    a = acc(t[0], q[0])
    for i in range(t.size - 1):
        dt = t[i + 1] - t[i]
        q[i + 1] = q[i] + v[i] * dt + 0.5 * a * dt**2
        a_new = acc(t[i + 1], q[i + 1])
        v[i + 1] = v[i] + 0.5 * (a + a_new) * dt
        a = a_new
    return q, v
