#pragma once

// Numerical flux functions at cell interfaces (Riemann solvers) and the
// analytic Euler flux.

#include <string>

#include "euler1d/eos.hpp"
#include "euler1d/types.hpp"

namespace euler1d {

Conserved physical_flux(const Primitive& W, const IdealGas& gas);

Conserved rusanov_flux(const Primitive& WL, const Primitive& WR, const IdealGas& gas);

// HLLC (Toro, Ch. 10) with Einfeldt/Roe-averaged wave speed estimates.
Conserved hllc_flux(const Primitive& WL, const Primitive& WR, const IdealGas& gas);

enum class FluxScheme { Rusanov, HLLC };

using FluxFunction = Conserved (*)(const Primitive&, const Primitive&, const IdealGas&);

FluxFunction select_flux(FluxScheme scheme);

// Case-insensitive; throws std::invalid_argument for unknown names.
FluxScheme parse_flux_scheme(const std::string& name);

}  // namespace euler1d
