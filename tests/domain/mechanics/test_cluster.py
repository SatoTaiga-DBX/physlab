"""ばねでつないだ質点: 力 = −∇U、内力の性質、対称性の破れの式。"""

import numpy as np
import pytest

from physlab.domain.mechanics.cluster import (
    FIELDS,
    ClusterParams,
    accel,
    angular_momentum,
    energy,
    external_energy,
    forces,
    momentum,
    predicted_rates,
    spring_energy,
)

Q = np.array([0.1, -0.3, 1.4, 0.2, 0.3, 1.1])
V = np.array([0.5, 0.2, -0.4, 0.7, 0.3, -0.9])
T = 0.8


def masses(p: ClusterParams) -> np.ndarray:
    return np.repeat(np.asarray(p.masses), 2)


@pytest.mark.parametrize("field", FIELDS)
def test_force_is_minus_gradient_of_potential(field: str) -> None:
    p = ClusterParams(field=field)

    def potential(q: np.ndarray) -> float:
        return float(spring_energy(q, p) + external_energy(T, q, p))

    h = 1e-6
    numeric = np.array(
        [-(potential(Q + h * e) - potential(Q - h * e)) / (2 * h) for e in np.eye(Q.size)]
    )
    assert masses(p) * accel(p)(T, Q) == pytest.approx(numeric, rel=1e-6, abs=1e-8)


def test_internal_forces_have_zero_sum_and_zero_torque() -> None:
    internal, external = forces(T, Q, ClusterParams())
    r = Q.reshape(-1, 2)
    assert internal.sum(axis=0) == pytest.approx(0.0, abs=1e-12)
    torque = np.sum(r[:, 0] * internal[:, 1] - r[:, 1] * internal[:, 0])
    assert torque == pytest.approx(0.0, abs=1e-12)
    assert external == pytest.approx(0.0)


@pytest.mark.parametrize("field", FIELDS)
def test_predicted_rates_match_the_time_derivative_along_the_motion(field: str) -> None:
    """dE/dt, dP/dt, dL/dt の式を、運動に沿った中心差分（独立な計算）と比べる。"""
    p = ClusterParams(field=field)
    a = accel(p)(T, Q)
    h = 1e-5

    def state(s: float) -> tuple[float, np.ndarray, np.ndarray]:
        return T + s, Q + s * V + 0.5 * s**2 * a, V + s * a

    def quantities(s: float) -> np.ndarray:
        t, q, v = state(s)
        return np.concatenate([[energy(t, q, v, p)], momentum(v, p), [angular_momentum(q, v, p)]])

    numeric = (quantities(h) - quantities(-h)) / (2 * h)
    d_e, d_p, d_l = predicted_rates(np.float64(T), Q, p)
    predicted = np.concatenate([[d_e], d_p, [d_l]])
    assert numeric == pytest.approx(predicted, abs=1e-6)


def test_unknown_field_is_rejected() -> None:
    with pytest.raises(ValueError, match="field"):
        ClusterParams(field="magnet")
