"""Exact solver for the 1D Euler Riemann problem with an ideal-gas equation of state.

Follows Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics", 3rd ed., Ch. 4.
States are primitive triples ``(rho, u, p)``.
"""
from __future__ import annotations

import numpy as np

State = tuple[float, float, float]

_TOL = 1e-12
_MAX_ITER = 100


def _f_and_df(p: float, state: State, gamma: float) -> tuple[float, float]:
    """Pressure function f_K(p) and its derivative for the side with ``state``."""
    rho_k, _, p_k = state
    a_k = np.sqrt(gamma * p_k / rho_k)
    if p > p_k:  # shock
        big_a = 2.0 / ((gamma + 1.0) * rho_k)
        big_b = (gamma - 1.0) / (gamma + 1.0) * p_k
        root = np.sqrt(big_a / (p + big_b))
        f = (p - p_k) * root
        df = root * (1.0 - (p - p_k) / (2.0 * (p + big_b)))
    else:  # rarefaction
        expo = (gamma - 1.0) / (2.0 * gamma)
        f = 2.0 * a_k / (gamma - 1.0) * ((p / p_k) ** expo - 1.0)
        df = (1.0 / (rho_k * a_k)) * (p / p_k) ** (-(gamma + 1.0) / (2.0 * gamma))
    return float(f), float(df)


def star_state(WL: State, WR: State, gamma: float) -> tuple[float, float]:
    """Return the star-region pressure and velocity ``(p_star, u_star)``.

    Raises ValueError if the data generate a vacuum and RuntimeError if the Newton
    iteration does not converge.
    """
    rho_l, u_l, p_l = WL
    rho_r, u_r, p_r = WR
    a_l = np.sqrt(gamma * p_l / rho_l)
    a_r = np.sqrt(gamma * p_r / rho_r)
    if 2.0 * (a_l + a_r) / (gamma - 1.0) <= u_r - u_l:
        raise ValueError("vacuum generated")

    p = max(_TOL, 0.5 * (p_l + p_r) - 0.125 * (u_r - u_l) * (rho_l + rho_r) * (a_l + a_r))
    for _ in range(_MAX_ITER):
        f_l, df_l = _f_and_df(p, WL, gamma)
        f_r, df_r = _f_and_df(p, WR, gamma)
        p_new = max(_TOL, p - (f_l + f_r + (u_r - u_l)) / (df_l + df_r))
        change = abs(p_new - p) / (0.5 * (p_new + p))
        p = p_new
        if change < _TOL:
            break
    else:
        raise RuntimeError("Newton iteration for the star pressure did not converge")

    f_l, _ = _f_and_df(p, WL, gamma)
    f_r, _ = _f_and_df(p, WR, gamma)
    return p, 0.5 * (u_l + u_r) + 0.5 * (f_r - f_l)


def solve(
    WL: State, WR: State, gamma: float, x: np.ndarray, t: float, x0: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sample the exact solution at positions ``x`` and time ``t > 0``.

    The initial discontinuity sits at ``x0``. Returns arrays ``(rho, u, p)``.
    """
    g = gamma
    rho_l, u_l, p_l = WL
    rho_r, u_r, p_r = WR
    a_l = np.sqrt(g * p_l / rho_l)
    a_r = np.sqrt(g * p_r / rho_r)
    p_s, u_s = star_state(WL, WR, g)

    s = (np.asarray(x, dtype=float) - x0) / t
    rho = np.empty_like(s)
    u = np.empty_like(s)
    p = np.empty_like(s)

    def assign(mask: np.ndarray, r, v, q) -> None:
        rho[mask], u[mask], p[mask] = r, v, q

    gm = (g - 1.0) / (g + 1.0)
    left = s <= u_s

    # Left of the contact.
    if p_s <= p_l:  # left rarefaction
        s_head = u_l - a_l
        a_star = a_l * (p_s / p_l) ** ((g - 1.0) / (2.0 * g))
        s_tail = u_s - a_star
        assign(left & (s <= s_head), rho_l, u_l, p_l)
        assign(left & (s > s_tail), rho_l * (p_s / p_l) ** (1.0 / g), u_s, p_s)
        fan = left & (s > s_head) & (s <= s_tail)
        sf = s[fan]
        c = 2.0 / (g + 1.0) * (a_l + 0.5 * (g - 1.0) * (u_l - sf))
        assign(fan, rho_l * (c / a_l) ** (2.0 / (g - 1.0)),
               2.0 / (g + 1.0) * (a_l + 0.5 * (g - 1.0) * u_l + sf),
               p_l * (c / a_l) ** (2.0 * g / (g - 1.0)))
    else:  # left shock
        pm = p_s / p_l
        s_shock = u_l - a_l * np.sqrt((g + 1.0) / (2.0 * g) * pm + (g - 1.0) / (2.0 * g))
        assign(left & (s <= s_shock), rho_l, u_l, p_l)
        assign(left & (s > s_shock), rho_l * (pm + gm) / (pm * gm + 1.0), u_s, p_s)

    # Right of the contact.
    right = ~left
    if p_s > p_r:  # right shock
        pm = p_s / p_r
        s_shock = u_r + a_r * np.sqrt((g + 1.0) / (2.0 * g) * pm + (g - 1.0) / (2.0 * g))
        assign(right & (s >= s_shock), rho_r, u_r, p_r)
        assign(right & (s < s_shock), rho_r * (pm + gm) / (pm * gm + 1.0), u_s, p_s)
    else:  # right rarefaction
        s_head = u_r + a_r
        a_star = a_r * (p_s / p_r) ** ((g - 1.0) / (2.0 * g))
        s_tail = u_s + a_star
        assign(right & (s >= s_head), rho_r, u_r, p_r)
        assign(right & (s <= s_tail), rho_r * (p_s / p_r) ** (1.0 / g), u_s, p_s)
        fan = right & (s < s_head) & (s > s_tail)
        sf = s[fan]
        c = 2.0 / (g + 1.0) * (a_r - 0.5 * (g - 1.0) * (u_r - sf))
        assign(fan, rho_r * (c / a_r) ** (2.0 / (g - 1.0)),
               2.0 / (g + 1.0) * (-a_r + 0.5 * (g - 1.0) * u_r + sf),
               p_r * (c / a_r) ** (2.0 * g / (g - 1.0)))

    return rho, u, p
