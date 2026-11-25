"""対称性と保存量: 壊した対称性に対応する量だけが、予測どおりに変化する。"""

import numpy as np
import pytest

from physlab.application.mechanics.symmetry import (
    center_of_mass_invariant,
    compare_integrators,
    rotating_rod_run,
    symmetry_experiment,
)

# 外場ごとに保存するはずの量（教科書の表と同じ）
CONSERVED = {
    "none": {"E", "Px", "Py", "L"},
    "gravity": {"E", "Px"},
    "trap": {"E", "L"},
    "pulsing_trap": {"L"},
}


@pytest.fixture(scope="module")
def runs():
    return symmetry_experiment()


def series(run) -> dict[str, np.ndarray]:
    return {"E": run.energy, "Px": run.momentum[:, 0], "Py": run.momentum[:, 1], "L": run.angular}


@pytest.mark.parametrize("field", list(CONSERVED))
def test_only_the_quantities_with_a_symmetry_are_conserved(runs, field: str) -> None:
    for name, value in series(runs[field]).items():
        change = np.max(np.abs(value - value[0]))
        if name in CONSERVED[field]:
            # E は速度ベルレの O(h²) の振動、運動量と角運動量は丸め誤差
            assert change < (2e-5 if name == "E" else 1e-10), (field, name, change)
        else:
            assert change > 0.5, (field, name, change)


@pytest.mark.parametrize("field", list(CONSERVED))
def test_changes_agree_with_the_rates_predicted_from_symmetry_breaking(runs, field: str) -> None:
    r = runs[field]
    assert r.energy == pytest.approx(r.predicted_energy, abs=2e-5)
    assert r.momentum == pytest.approx(r.predicted_momentum, abs=1e-8)
    assert r.angular == pytest.approx(r.predicted_angular, abs=1e-8)


def test_gravity_changes_vertical_momentum_linearly(runs) -> None:
    r = runs["gravity"]
    mass = r.params.total_mass
    assert r.momentum[:, 1] == pytest.approx(r.momentum[0, 1] - mass * r.params.g * r.t, abs=1e-8)


def test_free_system_keeps_m_r_minus_p_t_constant(runs) -> None:
    """ガリレイ変換の対称性に対応する保存量（重心の定理）。"""
    n = center_of_mass_invariant(runs["none"])
    assert n == pytest.approx(np.tile(n[0], (n.shape[0], 1)), abs=1e-10)


def test_verlet_keeps_momenta_to_roundoff_but_rk4_drifts_in_angular_momentum() -> None:
    c = compare_integrators()
    for method in ("verlet", "rk4"):
        assert c.momentum_error[method].max() < 1e-11  # 全運動量は 1 次の量: どちらも保存
    assert c.angular_error["verlet"].max() < 1e-11
    assert c.angular_error["rk4"].max() > 1e-9  # 2 次の量は RK4 では保存しない
    assert c.energy_error["verlet"].max() < 1e-3  # 有界に振動
    assert c.energy_error["rk4"].max() < 1e-5  # 小さいが一方向にずれる
    assert c.energy_error["rk4"][-1] > 0.5 * c.energy_error["rk4"].max()


def test_rotating_rod_conserves_energy_function_while_kinetic_energy_grows() -> None:
    rod = rotating_rod_run()
    assert rod.r == pytest.approx(rod.exact_r, rel=1e-9)
    assert np.max(np.abs(rod.energy_function - rod.energy_function[0])) < 1e-12
    assert rod.kinetic[-1] > 5 * rod.kinetic[0]
    # T の増加はモーターがした仕事に等しい
    assert rod.kinetic - rod.kinetic[0] == pytest.approx(rod.motor_work, abs=1e-5)
