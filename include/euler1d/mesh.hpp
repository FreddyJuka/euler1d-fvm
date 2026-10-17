#pragma once

// Uniform 1D cell-centered mesh with ghost cells.

namespace euler1d {

class Grid {
public:
    // Throws std::invalid_argument if n_cells < 1, x_max <= x_min, or n_ghost < 1.
    Grid(int n_cells, double x_min, double x_max, int n_ghost = 2);

    int n_cells() const { return n_cells_; }
    int n_ghost() const { return n_ghost_; }
    int size() const { return n_cells_ + 2 * n_ghost_; }
    int begin() const { return n_ghost_; }
    int end() const { return n_ghost_ + n_cells_; }
    double dx() const { return dx_; }
    double x_min() const { return x_min_; }
    double x_max() const { return x_max_; }

    // Cell center for storage index i (ghost cells included).
    double x_center(int i) const { return x_min_ + (i - n_ghost_ + 0.5) * dx_; }

private:
    int n_cells_;
    int n_ghost_;
    double x_min_;
    double x_max_;
    double dx_;
};

}  // namespace euler1d
