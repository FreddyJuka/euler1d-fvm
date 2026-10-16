#include "euler1d/flux.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace euler1d {

Conserved physical_flux(const Primitive& W, const IdealGas& gas) {
    const Conserved U = gas.to_conserved(W);
    return {U.mom, U.mom * W.u + W.p, W.u * (U.E + W.p)};
}

Conserved rusanov_flux(const Primitive& WL, const Primitive& WR, const IdealGas& gas) {
    const double s_max = std::max(std::abs(WL.u) + gas.sound_speed(WL),
                                  std::abs(WR.u) + gas.sound_speed(WR));
    const Conserved UL = gas.to_conserved(WL);
    const Conserved UR = gas.to_conserved(WR);
    return 0.5 * (physical_flux(WL, gas) + physical_flux(WR, gas)) - (0.5 * s_max) * (UR - UL);
}

namespace {

// Star state U*K for wave speed SK and contact speed S_star.
Conserved hllc_star_state(const Primitive& W, const Conserved& U, double SK, double S_star) {
    const double c = W.rho * (SK - W.u) / (SK - S_star);
    return {c, c * S_star,
            c * (U.E / W.rho + (S_star - W.u) * (S_star + W.p / (W.rho * (SK - W.u))))};
}

}  // namespace

Conserved hllc_flux(const Primitive& WL, const Primitive& WR, const IdealGas& gas) {
    const Conserved UL = gas.to_conserved(WL);
    const Conserved UR = gas.to_conserved(WR);

    const double HL = (UL.E + WL.p) / WL.rho;
    const double HR = (UR.E + WR.p) / WR.rho;
    const double sL = std::sqrt(WL.rho);
    const double sR = std::sqrt(WR.rho);
    const double u_t = (sL * WL.u + sR * WR.u) / (sL + sR);
    const double H_t = (sL * HL + sR * HR) / (sL + sR);
    const double a_t = std::sqrt((gas.gamma() - 1.0) * (H_t - 0.5 * u_t * u_t));

    const double SL = std::min(WL.u - gas.sound_speed(WL), u_t - a_t);
    const double SR = std::max(WR.u + gas.sound_speed(WR), u_t + a_t);

    if (SL >= 0.0) return physical_flux(WL, gas);
    if (SR <= 0.0) return physical_flux(WR, gas);

    const double dL = WL.rho * (SL - WL.u);
    const double dR = WR.rho * (SR - WR.u);
    const double S_star = (WR.p - WL.p + dL * WL.u - dR * WR.u) / (dL - dR);

    if (S_star >= 0.0) {
        const Conserved UsL = hllc_star_state(WL, UL, SL, S_star);
        return physical_flux(WL, gas) + SL * (UsL - UL);
    }
    const Conserved UsR = hllc_star_state(WR, UR, SR, S_star);
    return physical_flux(WR, gas) + SR * (UsR - UR);
}

FluxFunction select_flux(FluxScheme scheme) {
    switch (scheme) {
        case FluxScheme::Rusanov:
            return &rusanov_flux;
        case FluxScheme::HLLC:
            return &hllc_flux;
    }
    throw std::invalid_argument("select_flux: unknown FluxScheme");
}

FluxScheme parse_flux_scheme(const std::string& name) {
    std::string lower(name);
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lower == "rusanov") return FluxScheme::Rusanov;
    if (lower == "hllc") return FluxScheme::HLLC;
    throw std::invalid_argument("Unknown flux scheme '" + name +
                                "'. Valid names: rusanov, hllc");
}

}  // namespace euler1d
