#pragma once

// Equation of state for an ideal gas: conversions between conservative and
// primitive variables, pressure, sound speed, and total energy.

#include "euler1d/types.hpp"

namespace euler1d {

class IdealGas {
public:
    // Throws std::invalid_argument if gamma <= 1.
    explicit IdealGas(double gamma);

    double gamma() const { return gamma_; }

    // Neither conversion checks admissibility; use is_admissible().
    Primitive to_primitive(const Conserved& U) const;
    Conserved to_conserved(const Primitive& W) const;

    double sound_speed(const Primitive& W) const;

    // True iff rho > 0 and p > 0.
    bool is_admissible(const Conserved& U) const;

private:
    double gamma_;
};

}  // namespace euler1d
