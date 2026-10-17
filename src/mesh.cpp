#include "euler1d/mesh.hpp"

#include <stdexcept>

namespace euler1d {

Grid::Grid(int n_cells, double x_min, double x_max, int n_ghost)
    : n_cells_(n_cells), n_ghost_(n_ghost), x_min_(x_min), x_max_(x_max), dx_(0.0) {
    if (n_cells < 1) throw std::invalid_argument("Grid: n_cells must be >= 1");
    if (!(x_max > x_min)) throw std::invalid_argument("Grid: x_max must be > x_min");
    if (n_ghost < 1) throw std::invalid_argument("Grid: n_ghost must be >= 1");
    dx_ = (x_max - x_min) / n_cells;
}

}  // namespace euler1d
