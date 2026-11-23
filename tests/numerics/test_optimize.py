"""最急降下法: 2 次関数の最小点への収束と、途中経過の記録。"""

import numpy as np
import pytest

from physlab.numerics.optimize import steepest_descent

A = np.array([[3.0, 1.0], [1.0, 2.0]])
B = np.array([1.0, -2.0])


def grad(x: np.ndarray) -> np.ndarray:
    return A @ x - B  # f(x) = ½ x·A x − B·x の勾配


def test_converges_to_the_minimum_of_a_quadratic():
    x, _ = steepest_descent(grad, np.array([5.0, 5.0]), step=0.2, n_iter=300)
    assert x == pytest.approx(np.linalg.solve(A, B), abs=1e-10)


def test_records_requested_iterations_and_copies_them():
    x0 = np.array([5.0, 5.0])
    x, snaps = steepest_descent(grad, x0, step=0.2, n_iter=10, record_at=[0, 3, 10])
    assert sorted(snaps) == [0, 3, 10]
    assert snaps[0] == pytest.approx(x0)
    assert snaps[10] == pytest.approx(x)
    assert x0 == pytest.approx([5.0, 5.0])  # 入力は書き換えない
    assert not np.shares_memory(snaps[3], x)


def test_too_large_step_diverges():
    x, _ = steepest_descent(grad, np.array([1.0, 1.0]), step=1.0, n_iter=50)  # 1 > 2/λmax
    assert np.linalg.norm(x) > 1e3
