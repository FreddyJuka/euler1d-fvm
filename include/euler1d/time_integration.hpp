#pragma once

// Time advancement schemes and CFL-based time step selection.

#include <vector>

#include "euler1d/eos.hpp"
#include "euler1d/mesh.hpp"
#include "euler1d/types.hpp"

namespace euler1d {

// dt = cfl*dx / max over interior cells of (|u| + a).
// Throws std::invalid_argument if cfl <= 0 or cfl > 1.
double compute_dt(const std::vector<Primitive>& W, const Grid& grid, const IdealGas& gas,
                  double cfl);

}  // namespace euler1d
