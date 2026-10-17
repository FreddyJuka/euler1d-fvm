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


# ---- Sampling tests -------------------------------------------------------------------

TORO_TESTS = {
    1: ((1.0, 0.0, 1.0), (0.125, 0.0, 0.1), 0.42632, 0.26557),
    2: ((1.0, -2.0, 0.4), (1.0, 2.0, 0.4), 0.02185, 0.02185),
    3: ((1.0, 0.0, 1000.0), (1.0, 0.0, 0.01), 0.57506, 5.99924),
    4: ((1.0, 0.0, 0.01), (1.0, 0.0, 100.0), 5.99242, 0.57511),
    5: ((5.99924, 19.5975, 460.894), (5.99242, -6.19633, 46.0950), 14.2823, 31.0426),
}


@pytest.mark.parametrize("test_id", sorted(TORO_TESTS))
def test_star_densities_table_4_3(test_id):
    WL, WR, rho_l_ref, rho_r_ref = TORO_TESTS[test_id]
    t, x0 = 0.1, 0.5
    _, u_s = star_state(WL, WR, G)
    xc = x0 + u_s * t
    rho, _, _ = solve(WL, WR, G, np.array([xc - 1e-6, xc + 1e-6]), t, x0)
    assert rho[0] == pytest.approx(rho_l_ref, rel=1e-4)
    assert rho[1] == pytest.approx(rho_r_ref, rel=1e-4)


SOD = ((1.0, 0.0, 1.0), (0.125, 0.0, 0.1))


def _flux(rho, u, p):
    E = p / (G - 1.0) + 0.5 * rho * u * u
    return np.array([rho * u, rho * u * u + p, u * (E + p)]), np.array([rho, rho * u, E])


def test_sod_rankine_hugoniot_right_shock():
    WL, WR = SOD
    t, x0 = 0.2, 0.5
    _, u_s = star_state(WL, WR, G)
    x = np.linspace(x0 + (u_s + 0.1) * t, 1.0, 500001)  # right of the contact
    rho, u, p = solve(WL, WR, G, x, t, x0)
    k = int(np.argmax(np.abs(np.diff(rho))))
    x_shock = 0.5 * (x[k] + x[k + 1])
    speed = (x_shock - x0) / t
    # Cells k and k+1 sit just on either side of the jump (spacing 1e-6).
    F_l, U_l = _flux(rho[k], u[k], p[k])
    F_r, U_r = _flux(rho[k + 1], u[k + 1], p[k + 1])
    for i in range(3):
        assert speed * (U_r[i] - U_l[i]) == pytest.approx(F_r[i] - F_l[i], rel=1e-3)


def test_sod_fan_continuity():
    WL, WR = SOD
    t, x0 = 0.2, 0.5
    p_s, u_s = star_state(WL, WR, G)
    a_l = np.sqrt(G * WL[2] / WL[0])
    a_star = a_l * (p_s / WL[2]) ** ((G - 1.0) / (2.0 * G))
    for xe in (x0 + (WL[1] - a_l) * t, x0 + (u_s - a_star) * t):
        lo = np.array(solve(WL, WR, G, np.array([xe - 1e-9]), t, x0))
        hi = np.array(solve(WL, WR, G, np.array([xe + 1e-9]), t, x0))
        assert np.max(np.abs(hi - lo)) < 1e-6


def test_sod_fan_self_similarity():
    WL, WR = SOD
    t, x0 = 0.1, 0.5
    p_s, u_s = star_state(WL, WR, G)
    a_l = np.sqrt(G * WL[2] / WL[0])
    a_star = a_l * (p_s / WL[2]) ** ((G - 1.0) / (2.0 * G))
    x = x0 + np.linspace((WL[1] - a_l) * t, (u_s - a_star) * t, 20)[1:-1]
    a = np.array(solve(WL, WR, G, x, t, x0))
    b = np.array(solve(WL, WR, G, x0 + 2.0 * (x - x0), 2.0 * t, x0))
    assert np.max(np.abs(a - b)) < 1e-12
