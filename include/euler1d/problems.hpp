#pragma once

// Initial conditions for the supported test problems and exact-solution error norms.

#include <string>
#include <vector>

#include "euler1d/eos.hpp"
#include "euler1d/mesh.hpp"
#include "euler1d/types.hpp"

namespace euler1d {

enum class ProblemType { Riemann, DensityWave };

// Case-insensitive ("riemann", "density_wave"); throws std::invalid_argument otherwise.
ProblemType parse_problem_type(const std::string& name);

struct RiemannParams {
    double x0;
    double rho_L, u_L, p_L;
    double rho_R, u_R, p_R;
};

struct DensityWaveParams {
    double rho0;
    double amplitude;
    double u0;
    double p0;
    int wavenumber;
};

struct ProblemConfig {
    ProblemType type;
    RiemannParams riemann;
    DensityWaveParams density_wave;
};

struct ErrorNorms {
    double L1;
    double L2;
    double Linf;
};

// Exact cell average of the advected density wave over the cell centered at xc.
double density_wave_cell_average(const DensityWaveParams& params, const Grid& grid, double xc,
                                 double t);

// Fills interior cells (result has size grid.size(); ghost cells are zero).
// Throws std::invalid_argument for invalid density-wave parameters.
std::vector<Conserved> make_initial_condition(const ProblemConfig& problem, const Grid& grid,
                                              const IdealGas& gas);

// Density error against the exact cell averages over interior cells.
ErrorNorms density_wave_errors(const std::vector<Conserved>& U, const Grid& grid,
                               const DensityWaveParams& params, double t);

}  // namespace euler1d
