import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))

from exact_riemann import solve, star_state  # noqa: E402

G = 1.4

# (WL, WR, p*, u*, rtol for p*) -- Toro 3rd ed., Table 4.3. Test 2 lists p* with only three
# significant digits (0.00189), so its tolerance reflects the tabulated precision.
TABLE_4_3 = [
    ((1.0, 0.0, 1.0), (0.125, 0.0, 0.1), 0.30313, 0.92745, 1e-4),
    ((1.0, -2.0, 0.4), (1.0, 2.0, 0.4), 0.00189, 0.0, 5e-3),
    ((1.0, 0.0, 1000.0), (1.0, 0.0, 0.01), 460.894, 19.5975, 1e-4),
    ((1.0, 0.0, 0.01), (1.0, 0.0, 100.0), 46.0950, -6.19633, 1e-4),
    ((5.99924, 19.5975, 460.894), (5.99242, -6.19633, 46.0950), 1691.64, 8.68975, 1e-4),
]


@pytest.mark.parametrize("WL,WR,p_ref,u_ref,rtol", TABLE_4_3)
def test_star_state_table_4_3(WL, WR, p_ref, u_ref, rtol):
    p, u = star_state(WL, WR, G)
    assert p == pytest.approx(p_ref, rel=rtol)
    if u_ref == 0.0:
        assert abs(u) < 1e-10
    else:
        assert u == pytest.approx(u_ref, rel=1e-4)


@pytest.mark.parametrize("WL,WR,p_ref,u_ref,rtol", TABLE_4_3)
def test_far_field(WL, WR, p_ref, u_ref, rtol):
    x = np.array([-1e6, 1e6])
    rho, u, p = solve(WL, WR, G, x, 0.1, 0.0)
    assert (rho[0], u[0], p[0]) == WL
    assert (rho[1], u[1], p[1]) == WR


def test_vacuum_raises():
    with pytest.raises(ValueError, match="vacuum generated"):
        star_state((1.0, -10.0, 1.0), (1.0, 10.0, 1.0), G)


def test_uniform_state():
    W = (0.7, 0.3, 1.2)
    x = np.linspace(0.0, 1.0, 101)
    rho, u, p = solve(W, W, G, x, 0.2, 0.5)
    np.testing.assert_allclose(rho, W[0], rtol=1e-13)
    np.testing.assert_allclose(u, W[1], rtol=1e-13)
    np.testing.assert_allclose(p, W[2], rtol=1e-13)
