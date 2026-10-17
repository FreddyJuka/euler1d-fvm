#!/usr/bin/env python3
"""Independent 50-digit HLLC reference fluxes (gamma = 1.4) for tests/test_flux.cpp.

Re-implements the HLLC flux with Einfeldt/Roe-averaged wave speeds from the formulas in
README.md using mpmath; deliberately imports no project code.
"""
from mpmath import mp, mpf, sqrt

mp.dps = 50
GAMMA = mpf("1.4")

CASES = {
    "Sod": ((1, 0, 1), (mpf("0.125"), 0, mpf("0.1"))),
    "GeneralSubsonic": ((mpf("1.0"), mpf("0.3"), mpf("1.0")), (mpf("0.5"), mpf("-0.2"), mpf("0.4"))),
}


def conserved(W):
    rho, u, p = (mpf(v) for v in W)
    return [rho, rho * u, p / (GAMMA - 1) + rho * u * u / 2]


def physical_flux(W):
    rho, u, p = (mpf(v) for v in W)
    U = conserved(W)
    return [rho * u, rho * u * u + p, u * (U[2] + p)]


def hllc_flux(WL, WR):
    rho_l, u_l, p_l = (mpf(v) for v in WL)
    rho_r, u_r, p_r = (mpf(v) for v in WR)
    UL, UR = conserved(WL), conserved(WR)
    a_l, a_r = sqrt(GAMMA * p_l / rho_l), sqrt(GAMMA * p_r / rho_r)
    H_l, H_r = (UL[2] + p_l) / rho_l, (UR[2] + p_r) / rho_r
    w_l, w_r = sqrt(rho_l), sqrt(rho_r)
    u_t = (w_l * u_l + w_r * u_r) / (w_l + w_r)
    H_t = (w_l * H_l + w_r * H_r) / (w_l + w_r)
    a_t = sqrt((GAMMA - 1) * (H_t - u_t * u_t / 2))
    S_l = min(u_l - a_l, u_t - a_t)
    S_r = max(u_r + a_r, u_t + a_t)
    S_s = (p_r - p_l + rho_l * u_l * (S_l - u_l) - rho_r * u_r * (S_r - u_r)) / (
        rho_l * (S_l - u_l) - rho_r * (S_r - u_r)
    )

    def star(rho, u, p, U, S):
        f = rho * (S - u) / (S - S_s)
        return [f, f * S_s, f * (U[2] / rho + (S_s - u) * (S_s + p / (rho * (S - u))))]

    FL, FR = physical_flux(WL), physical_flux(WR)
    if S_l >= 0:
        return FL
    if S_s >= 0:
        Us = star(rho_l, u_l, p_l, UL, S_l)
        return [FL[i] + S_l * (Us[i] - UL[i]) for i in range(3)]
    if S_r > 0:
        Us = star(rho_r, u_r, p_r, UR, S_r)
        return [FR[i] + S_r * (Us[i] - UR[i]) for i in range(3)]
    return FR


if __name__ == "__main__":
    for name, (WL, WR) in CASES.items():
        F = hllc_flux(WL, WR)
        print(f"{name}: " + ", ".join(mp.nstr(f, 16) for f in F))
