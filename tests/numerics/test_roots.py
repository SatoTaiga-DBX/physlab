"""2 分法。"""

import numpy as np
import pytest

from physlab.numerics.roots import bisect


def test_finds_roots_of_many_functions_at_once() -> None:
    a = np.array([1.0, 2.0, 9.0])
    root = bisect(lambda x: a - x**2, np.zeros(3), np.full(3, 10.0))
    assert root == pytest.approx(np.sqrt(a), rel=1e-14)
    assert np.all(a - root**2 >= -1e-15)  # 正の側の端を返す


def test_requires_a_bracket() -> None:
    with pytest.raises(ValueError):
        bisect(lambda x: 1.0 + x, np.zeros(1), np.ones(1))
