"""力学的相似の実験: 冪ポテンシャルの周期、軌道の重ね書き、ビリアル定理、ケプラーの法則。

設計書: docs/design/006-conservation-laws.md ／ 教科書: 第 4 章 §力学的相似
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid

from physlab.application.mechanics.pendulum import Trajectory, measure_period
from physlab.domain.mechanics import two_body
from physlab.domain.mechanics.power_law import accel, period_exact, potential, time_scale
from physlab.numerics.integrators import Array, velocity_verlet

KS = (1.0, 2.0, 4.0, 6.0)


def oscillate(
    k: float, amplitude: float, n_periods: float = 1.2, steps_per_period: int = 4000
) -> tuple[Array, Array, Array]:
    """``x(0) = A, ẋ(0) = 0`` を速度ベルレ法で解き ``(t, x, v)`` を返す。刻みは周期に比例させる。"""
    dt = period_exact(k, amplitude) / steps_per_period
    t = np.arange(int(round(n_periods * steps_per_period)) + 1) * dt
    x, v = velocity_verlet(accel(k), [amplitude], [0.0], t)
    return t, x[:, 0], v[:, 0]


# ---- 周期と振幅 --------------------------------------------------------------------------


@dataclass(frozen=True)
class PeriodSweep:
    amplitudes: Array
    measured: dict[float, Array]  # k -> 数値解の周期
    exact: dict[float, Array]  # k -> 公式
    slopes: dict[float, float]  # log T と log A の傾き（理論は 1 − k/2）


def sweep_periods(ks: tuple[float, ...] = KS, amplitudes: Array | None = None) -> PeriodSweep:
    amplitudes = np.geomspace(0.25, 4.0, 7) if amplitudes is None else np.asarray(amplitudes)
    measured, exact, slopes = {}, {}, {}
    for k in ks:
        periods = []
        for a in amplitudes:
            t, x, v = oscillate(k, a)
            periods.append(measure_period(Trajectory(t, x, v)))
        measured[k] = np.array(periods)
        exact[k] = np.array([period_exact(k, a) for a in amplitudes])
        slopes[k] = float(np.polyfit(np.log(amplitudes), np.log(measured[k]), 1)[0])
    return PeriodSweep(amplitudes, measured, exact, slopes)


@dataclass(frozen=True)
class KeplerSweep:
    """万有引力（U ∝ 1/r、k = −1）。長半径 a を大きさの目安にすると T ∝ a^{3/2}。"""

    semi_major: Array
    measured: Array
    exact: Array
    slope: float


def kepler_periods(periapses: tuple[float, ...] = (0.5, 0.8, 1.2, 2.0)) -> KeplerSweep:
    """重心系の二体運動で、近点から次の近点までの時間を測る。"""
    axes, measured, exact = [], [], []
    for periapsis in periapses:
        params = two_body.TwoBodyParams(periapsis=periapsis)
        dt = params.period / 20000
        t = np.arange(int(1.2 * 20000) + 1) * dt
        q, v = velocity_verlet(two_body.accel(params), *two_body.initial_state_cm(params), t)
        r = q[:, 0:2] - q[:, 2:4]
        u = v[:, 0:2] - v[:, 2:4]
        radial = np.sum(r * u, axis=1)  # r・v。近点で − から + に変わる
        i = np.nonzero((radial[1:-1] < 0) & (radial[2:] >= 0))[0][0] + 1
        measured.append(t[i] - radial[i] * dt / (radial[i + 1] - radial[i]))
        exact.append(params.period)
        axes.append(params.semi_major_axis)
    axes, measured = np.array(axes), np.array(measured)
    slope = float(np.polyfit(np.log(axes), np.log(measured), 1)[0])
    return KeplerSweep(axes, measured, np.array(exact), slope)


# ---- 相似な軌道の重ね書き -----------------------------------------------------------------


@dataclass(frozen=True)
class Collapse:
    k: float
    amplitudes: tuple[float, ...]
    t: list[Array]  # そのままの時刻
    x: list[Array]  # そのままの位置
    scaled_t: list[Array]  # t / A^{1−k/2}
    scaled_x: list[Array]  # x / A
    max_difference: float  # 振幅 1 の運動との差（相似変換した後）の最大


def collapse(
    k: float = 4.0, amplitudes: tuple[float, ...] = (0.5, 1.0, 2.0), n_periods: float = 2.0
) -> Collapse:
    """振幅の違う運動を相似変換すると、振幅 1 の運動に重なることを確かめる。"""
    reference_t, reference_x, _ = oscillate(k, 1.0, n_periods)
    ts, xs, st, sx, worst = [], [], [], [], 0.0
    for a in amplitudes:
        t, x, _ = oscillate(k, a, n_periods)
        scaled_t, scaled_x = t / time_scale(k, a), x / a
        worst = max(
            worst, float(np.max(np.abs(scaled_x - np.interp(scaled_t, reference_t, reference_x))))
        )
        ts.append(t)
        xs.append(x)
        st.append(scaled_t)
        sx.append(scaled_x)
    return Collapse(k, tuple(amplitudes), ts, xs, st, sx, worst)


# ---- ビリアル定理 ------------------------------------------------------------------------


@dataclass(frozen=True)
class Virial:
    ks: tuple[float, ...]
    periods: dict[float, Array]  # k -> 時刻 / 周期
    running: dict[float, Array]  # k -> 時刻 0 から t までの 2⟨T⟩/⟨U⟩
    over_whole_periods: dict[float, float]  # k -> 整数周期ぶんの 2⟨T⟩/⟨U⟩（理論は k）


def _upward_zero_times(t: Array, x: Array) -> Array:
    """``x`` が負から正に変わる時刻（線形補間）。隣り合う 2 つの間隔が数値解の 1 周期。"""
    i = np.nonzero((x[:-1] < 0) & (x[1:] >= 0))[0]
    return t[i] - x[i] * (t[i + 1] - t[i]) / (x[i + 1] - x[i])


def virial(ks: tuple[float, ...] = KS, n_periods: int = 8, steps_per_period: int = 4000) -> Virial:
    """時間平均の比 ``2⟨T⟩/⟨U⟩``。整数周期ぶんは、数値解自身の周期（ゼロ交差）で区切る。

    厳密な周期の整数倍で区切ると、力が不連続な ``k = 1`` では数値解の周期のわずかなずれ（約 0.1%）が
    半端な区間として残り、比が 1% ずれる。
    """
    times, running, whole = {}, {}, {}
    for k in ks:
        t, x, v = oscillate(k, 1.0, n_periods, steps_per_period)
        kinetic = cumulative_trapezoid(0.5 * v**2, t, initial=0.0)
        pot = cumulative_trapezoid(potential(x, k), t, initial=0.0)
        times[k] = t[1:] / period_exact(k, 1.0)
        running[k] = 2.0 * kinetic[1:] / pot[1:]
        zeros = _upward_zero_times(t, x)
        first, last = zeros[0], zeros[-1]
        d_kinetic = np.interp(last, t, kinetic) - np.interp(first, t, kinetic)
        d_pot = np.interp(last, t, pot) - np.interp(first, t, pot)
        whole[k] = float(2.0 * d_kinetic / d_pot)
    return Virial(tuple(ks), times, running, whole)


def kepler_virial(periapsis: float = 1.0, steps: int = 20000) -> float:
    """ケプラー運動の 1 周期ぶんの ``2⟨T⟩/⟨U⟩``（理論は k = −1）。"""
    params = two_body.TwoBodyParams(periapsis=periapsis)
    t = np.arange(steps + 1) * (params.period / steps)
    q, v = velocity_verlet(two_body.accel(params), *two_body.initial_state_cm(params), t)
    kinetic = np.trapezoid(two_body.kinetic_energy(v, params), t)
    pot = np.trapezoid(two_body.potential(q, params), t)
    return float(2.0 * kinetic / pot)
