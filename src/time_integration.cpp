#include "euler1d/time_integration.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace euler1d {

double compute_dt(const std::vector<Primitive>& W, const Grid& grid, const IdealGas& gas,
                  double cfl) {
    if (!(cfl > 0.0 && cfl <= 1.0)) {
        throw std::invalid_argument("compute_dt: cfl must be in (0, 1]");
    }
    double s_max = 0.0;
    for (int i = grid.begin(); i < grid.end(); ++i) {
        s_max = std::max(s_max, std::abs(W[i].u) + gas.sound_speed(W[i]));
    }
    return cfl * grid.dx() / s_max;
}

}  // namespace euler1d
