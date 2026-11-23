"""積分器の検証: 調和振動子 y'' = -y（解析解 y = cos t）で収束次数と保存則を確かめる。"""

import numpy as np
import pytest

from physlab.numerics.integrators import integrate, velocity_verlet


def oscillator(_t, y):
    return np.array([y[1], -y[0]])


def final_error(method: str, dt: float, t_end: float = 2.0) -> float:
    t = np.arange(0.0, t_end + dt / 2, dt)
    y = integrate(oscillator, [1.0, 0.0], t, method=method)
    return abs(y[-1, 0] - np.cos(t[-1]))


@pytest.mark.parametrize(("method", "order"), [("euler", 1), ("rk4", 4)])
def test_convergence_order(method: str, order: int) -> None:
    """刻みを半分にすると誤差は 2^order 分の 1 になる。"""
    ratio = final_error(method, 0.01) / final_error(method, 0.005)
    assert ratio == pytest.approx(2**order, rel=0.15)


def test_verlet_second_order_and_bounded_energy() -> None:
    def acc(_t, q):
        return -q

    def err(dt: float) -> float:
        t = np.arange(0.0, 2.0 + dt / 2, dt)
        q, _ = velocity_verlet(acc, [1.0], [0.0], t)
        return abs(q[-1, 0] - np.cos(t[-1]))

    assert err(0.01) / err(0.005) == pytest.approx(4, rel=0.15)

    # 長時間でもエネルギーがドリフトしない（振動はするが有界）
    t = np.arange(0.0, 1000.0, 0.05)
    q, v = velocity_verlet(acc, [1.0], [0.0], t)
    energy = 0.5 * v[:, 0] ** 2 + 0.5 * q[:, 0] ** 2
    assert np.abs(energy - 0.5).max() < 1e-3


def test_euler_energy_grows() -> None:
    """前進オイラーは調和振動子のエネルギーを単調に増やしてしまう。"""
    t = np.arange(0.0, 50.0, 0.01)
    y = integrate(oscillator, [1.0, 0.0], t, method="euler")
    energy = 0.5 * y[:, 1] ** 2 + 0.5 * y[:, 0] ** 2
    assert energy[-1] > 1.5 * energy[0]


def test_euler_energy_factor_per_step_is_exactly_1_plus_dt_squared() -> None:
    """調和振動子では 1 ステップごとにエネルギーがちょうど (1 + dt²) 倍になる（教科書 00-1 章）。"""
    dt = 0.1
    y0 = np.array([1.0, 0.0])
    y1 = integrate(oscillator, y0, [0.0, dt], method="euler")[-1]
    assert (y1 @ y1) / (y0 @ y0) == pytest.approx(1 + dt**2, rel=1e-12)
