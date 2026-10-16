#include "euler1d/eos.hpp"

#include <cmath>
#include <stdexcept>

namespace euler1d {

IdealGas::IdealGas(double gamma) : gamma_(gamma) {
    if (!(gamma > 1.0)) {
        throw std::invalid_argument("IdealGas: gamma must be > 1");
    }
}

Primitive IdealGas::to_primitive(const Conserved& U) const {
    const double u = U.mom / U.rho;
    const double p = (gamma_ - 1.0) * (U.E - 0.5 * U.rho * u * u);
    return {U.rho, u, p};
}

Conserved IdealGas::to_conserved(const Primitive& W) const {
    return {W.rho, W.rho * W.u, W.p / (gamma_ - 1.0) + 0.5 * W.rho * W.u * W.u};
}

double IdealGas::sound_speed(const Primitive& W) const {
    return std::sqrt(gamma_ * W.p / W.rho);
}

bool IdealGas::is_admissible(const Conserved& U) const {
    if (!(U.rho > 0.0)) return false;
    return to_primitive(U).p > 0.0;
}

}  // namespace euler1d
