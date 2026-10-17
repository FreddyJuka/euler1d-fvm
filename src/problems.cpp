#include "euler1d/problems.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace euler1d {

namespace {

const double kPi = 3.14159265358979323846;

void validate(const DensityWaveParams& p) {
    if (!(p.rho0 - std::abs(p.amplitude) > 0.0)) {
        throw std::invalid_argument("density_wave: rho0 - |amplitude| must be > 0");
    }
    if (!(p.p0 > 0.0)) throw std::invalid_argument("density_wave: p0 must be > 0");
    if (p.wavenumber < 1) throw std::invalid_argument("density_wave: wavenumber must be >= 1");
}

}  // namespace

ProblemType parse_problem_type(const std::string& name) {
    std::string lower(name);
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lower == "riemann") return ProblemType::Riemann;
    if (lower == "density_wave") return ProblemType::DensityWave;
    throw std::invalid_argument("Unknown problem '" + name +
                                "'. Valid names: riemann, density_wave");
}

double density_wave_cell_average(const DensityWaveParams& params, const Grid& grid, double xc,
                                 double t) {
    const double L = grid.x_max() - grid.x_min();
    const double kappa = 2.0 * kPi * params.wavenumber / L;
    const double half = 0.5 * kappa * grid.dx();
    return params.rho0 + params.amplitude * std::sin(kappa * (xc - params.u0 * t - grid.x_min())) *
                             (std::sin(half) / half);
}

std::vector<Conserved> make_initial_condition(const ProblemConfig& problem, const Grid& grid,
                                              const IdealGas& gas) {
    std::vector<Conserved> U(grid.size(), Conserved{0.0, 0.0, 0.0});
    if (problem.type == ProblemType::Riemann) {
        const RiemannParams& r = problem.riemann;
        const Conserved UL = gas.to_conserved({r.rho_L, r.u_L, r.p_L});
        const Conserved UR = gas.to_conserved({r.rho_R, r.u_R, r.p_R});
        for (int i = grid.begin(); i < grid.end(); ++i) {
            U[i] = grid.x_center(i) < r.x0 ? UL : UR;
        }
    } else {
        const DensityWaveParams& w = problem.density_wave;
        validate(w);
        for (int i = grid.begin(); i < grid.end(); ++i) {
            const double rho = density_wave_cell_average(w, grid, grid.x_center(i), 0.0);
            U[i] = {rho, w.u0 * rho, w.p0 / (gas.gamma() - 1.0) + 0.5 * w.u0 * w.u0 * rho};
        }
    }
    return U;
}

ErrorNorms density_wave_errors(const std::vector<Conserved>& U, const Grid& grid,
                               const DensityWaveParams& params, double t) {
    double l1 = 0.0, l2 = 0.0, linf = 0.0;
    for (int i = grid.begin(); i < grid.end(); ++i) {
        const double e =
            U[i].rho - density_wave_cell_average(params, grid, grid.x_center(i), t);
        l1 += std::abs(e) * grid.dx();
        l2 += e * e * grid.dx();
        linf = std::max(linf, std::abs(e));
    }
    return {l1, std::sqrt(l2), linf};
}

}  // namespace euler1d
