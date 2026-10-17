#pragma once

// Top-level solver driver: owns the grid and state and orchestrates
// reconstruction, flux evaluation, boundary conditions, and time integration.

#include <array>
#include <functional>
#include <vector>

#include "euler1d/boundary.hpp"
#include "euler1d/eos.hpp"
#include "euler1d/flux.hpp"
#include "euler1d/mesh.hpp"
#include "euler1d/reconstruction.hpp"
#include "euler1d/types.hpp"

namespace euler1d {

struct SolverConfig {
    FluxScheme flux;
    ReconstructionConfig recon;
    BoundaryType left;
    BoundaryType right;
    double cfl;
    double t_final;
    int max_steps;
};

struct StepInfo {
    int step;
    double t;
    double dt;
    // sqrt(sum_interior ((U^{n+1} - U^n)/dt)^2 * dx) for rho, mom, E.
    std::array<double, 3> residual_L2;
};

struct RunStats {
    int n_steps;
    double t_final_reached;
};

class Solver {
public:
    // U0.size() must equal grid.size() (ghost values are ignored). Throws
    // std::invalid_argument for invalid boundaries, cfl not in (0,1], t_final not finite or <= 0,
    // max_steps <= 0, or a size mismatch.
    Solver(const Grid& grid, const IdealGas& gas, const SolverConfig& cfg,
           std::vector<Conserved> U0);

    // Advances one SSP-RK3 step with dt = min(compute_dt, dt_max); returns dt used.
    // Throws std::invalid_argument if !(dt_max > 0), std::runtime_error if any interior cell
    // becomes inadmissible.
    double step(double dt_max);

    // Runs until t_final or max_steps; calls the observer after each step if set.
    RunStats run(const std::function<void(const StepInfo&)>& observer = {});

    double time() const { return t_; }
    const std::vector<Conserved>& state() const { return U_; }
    const Grid& grid() const { return grid_; }

private:
    // out_interior = L(U_in); U_in is copied, so its ghost cells are ignored.
    void spatial_operator(const std::vector<Conserved>& U_in, std::vector<Conserved>& L);
    void check_admissible() const;

    Grid grid_;
    IdealGas gas_;
    SolverConfig cfg_;
    FluxFunction flux_;

    std::vector<Conserved> U_;
    std::vector<Conserved> U1_, U2_, L_, Uw_;
    std::vector<Primitive> W_, WL_, WR_;
    std::vector<Conserved> F_;

    double t_ = 0.0;
    int step_count_ = 0;
    std::array<double, 3> last_residual_{};
};

}  // namespace euler1d
