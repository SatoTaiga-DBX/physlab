"""最小作用の原理: 作用が ε の 2 次式であること、最小か鞍点か、緩和が真の運動に収束すること。"""

import numpy as np
import pytest

from physlab.application.mechanics.least_action import (
    action_landscape,
    relax_path,
    scan_action,
    scan_curvature_vs_duration,
)
from physlab.domain.mechanics.lagrangian import (
    harmonic_oscillator,
    uniform_gravity,
)

GRAVITY = uniform_gravity(1.0, 9.81)
OSC = harmonic_oscillator(1.0, 2 * np.pi)  # π/ω = 0.5 s
EPS = np.linspace(-2.0, 2.0, 21)


@pytest.mark.parametrize(
    ("system", "tau", "q1", "q2", "mode"),
    [
        (GRAVITY, 2.0, 0.0, 0.0, 1),
        (GRAVITY, 2.0, 0.0, 3.0, 2),
        (OSC, 0.3, 1.0, 0.0, 1),
        (OSC, 0.8, 1.0, 0.0, 1),
        (OSC, 0.8, 1.0, 0.0, 2),
    ],
)
def test_action_is_exactly_quadratic_with_the_predicted_curvature(system, tau, q1, q2, mode):
    scan = scan_action(system, 0.0, q1, tau, q2, EPS, mode=mode)
    assert scan.fit_curvature == pytest.approx(scan.theory_curvature, rel=1e-5)
    # 停留: ε の 1 次の項は 0（2 次の項に比べて十分小さい）
    assert abs(scan.fit_linear) < 1e-5 * abs(scan.theory_curvature)
    # S(ε) = S* + c ε²（許容は S の変化幅に対する割合。q̇ の差分近似の誤差 O(Δt²) を含む）
    expected = scan.action_true + scan.theory_curvature * EPS**2
    assert scan.action == pytest.approx(expected, abs=1e-5 * np.ptp(scan.action))


def test_gravity_action_of_the_true_path_matches_closed_form():
    scan = scan_action(GRAVITY, 0.0, 0.0, 2.0, 0.0, EPS)
    assert scan.action_true == pytest.approx(-(9.81**2) * 2.0**3 / 24, rel=1e-6)
    assert scan.fit_curvature == pytest.approx(0.25 * 2.0 * (np.pi / 2.0) ** 2, rel=1e-5)


def test_true_path_minimizes_action_for_a_short_interval():
    scan = scan_action(OSC, 0.0, 1.0, 0.3, 0.0, EPS)
    others = np.delete(scan.action, len(EPS) // 2)
    assert np.all(others > scan.action_true)


def test_true_path_is_only_a_saddle_for_a_long_interval():
    """τ > π/ω では、モード 1 の方向に作用が減る（最小ではない）。モード 2 では増える。"""
    mode1 = scan_action(OSC, 0.0, 1.0, 0.8, 0.0, EPS, mode=1)
    mode2 = scan_action(OSC, 0.0, 1.0, 0.8, 0.0, EPS, mode=2)
    nonzero = EPS != 0
    assert np.all(mode1.action[nonzero] < mode1.action_true)
    assert np.all(mode2.action[nonzero] > mode2.action_true)


def test_curvature_numerics_match_theory_and_change_sign_at_n_pi_over_omega():
    taus = np.linspace(0.2, 1.6, 141)
    scan = scan_curvature_vs_duration(OSC, taus, modes=(1, 2, 3))
    for n in (1, 2, 3):
        # 許容は c_n の最大値に対する割合（q̇ の差分近似の誤差 O(Δt²) を含む）
        assert scan.numeric[n] == pytest.approx(
            scan.theory[n], abs=1e-5 * np.max(np.abs(scan.theory[n]))
        )
        # 符号が + から − に変わる区間の線形補間が nπ/ω = 0.5 n
        c = scan.numeric[n]
        i = np.nonzero((c[:-1] > 0) & (c[1:] <= 0))[0][0]
        root = taus[i] - c[i] * (taus[i + 1] - taus[i]) / (c[i + 1] - c[i])
        assert root == pytest.approx(0.5 * n, abs=1e-3)


def test_landscape_is_a_bowl_for_short_and_a_saddle_for_long_interval():
    short = action_landscape(OSC, 0.0, 1.0, 0.3, 0.0)
    assert short.delta_action.min() == pytest.approx(0.0, abs=1e-6)
    assert np.all(short.delta_action >= -1e-6)  # どの方向にも上がる

    long = action_landscape(OSC, 0.0, 1.0, 0.8, 0.0)
    assert long.delta_action.min() < -1.0  # 下がる方向がある
    assert long.delta_action.max() > 1.0  # 上がる方向もある（鞍点）
    mid = len(long.a) // 2
    assert long.delta_action[mid, mid] == pytest.approx(0.0, abs=1e-6)
    # a 軸（b = 0）に沿っては c_1 a²
    on_a_axis = long.delta_action[mid, :]
    assert on_a_axis == pytest.approx(long.curvature_a * long.a**2, rel=1e-3, abs=1e-5)


def test_relaxation_lowers_the_action_monotonically_and_reaches_the_parabola():
    r = relax_path(GRAVITY, 0.0, 0.0, 2.0, 0.0)
    assert r.iterations[0] == 0 and r.iterations[-1] == 8000
    assert np.all(np.diff(r.action) <= 1e-9)  # 作用は下がり続ける
    assert r.action[0] > r.action[-1]
    assert r.action[-1] == pytest.approx(r.action_true, rel=1e-9)
    assert r.max_error[0] > 1.0  # 初めはまったく違う形
    assert r.max_error[-1] < 1e-6  # 最後は放物線に節点上で一致
    # 作用の最小値は解析的な S* = −g²τ³/24 にほぼ等しい（離散化の誤差の範囲）
    assert r.action[-1] == pytest.approx(-(9.81**2) * 2.0**3 / 24, rel=1e-3)


def test_relaxation_for_oscillator_converges_to_the_analytic_path():
    r = relax_path(OSC, 0.0, 1.0, 0.3, 0.0)
    assert r.max_error[-1] < 1e-3


def test_relaxation_refuses_a_saddle():
    with pytest.raises(ValueError, match="鞍点"):
        relax_path(OSC, 0.0, 1.0, 0.8, 0.0)
