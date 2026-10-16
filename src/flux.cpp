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

FluxFunction select_flux(FluxScheme scheme) {
    switch (scheme) {
        case FluxScheme::Rusanov:
            return &rusanov_flux;
    }
    throw std::invalid_argument("select_flux: unknown FluxScheme");
}

FluxScheme parse_flux_scheme(const std::string& name) {
    std::string lower(name);
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lower == "rusanov") return FluxScheme::Rusanov;
    throw std::invalid_argument("Unknown flux scheme '" + name + "'. Valid names: rusanov");
}

}  // namespace euler1d
