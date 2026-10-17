#include <cmath>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "euler1d/solver.hpp"
#include "euler1d/time_integration.hpp"

using namespace euler1d;

namespace {

const IdealGas gas(1.4);
const double kPi = 3.14159265358979323846;

std::vector<Conserved> make_state(const Grid& grid, Primitive (*init)(double)) {
    std::vector<Conserved> U(grid.size());
    for (int i = 0; i < grid.size(); ++i) U[i] = gas.to_conserved(init(grid.x_center(i)));
    return U;
}

std::array<double, 3> integral(const Solver& s) {
    std::array<double, 3> sum{0.0, 0.0, 0.0};
    const Grid& g = s.grid();
    for (int j = g.begin(); j < g.end(); ++j) {
        sum[0] += s.state()[j].rho * g.dx();
        sum[1] += s.state()[j].mom * g.dx();
        sum[2] += s.state()[j].E * g.dx();
    }
    return sum;
}

Primitive sine_wave(double x) { return {1.0 + 0.2 * std::sin(2.0 * kPi * x), 1.0, 1.0}; }
Primitive sod(double x) {
    return x < 0.5 ? Primitive{1.0, 0.0, 1.0} : Primitive{0.125, 0.0, 0.1};
}

SolverConfig sod_config(int max_steps = 100000) {
    return {FluxScheme::HLLC,
            {ReconstructionScheme::MUSCL, Limiter::Minmod},
            BoundaryType::Transmissive,
            BoundaryType::Transmissive,
            0.8,
            0.2,
            max_steps};
}

}  // namespace

TEST(TimeStep, ComputeDt) {
    const Grid grid(10, 0.0, 1.0);
    const std::vector<Primitive> W(grid.size(), Primitive{1.0, 0.5, 1.0});
    EXPECT_NEAR(compute_dt(W, grid, gas, 0.5), 0.5 * 0.1 / (0.5 + std::sqrt(1.4)), 1e-14);
    EXPECT_THROW(compute_dt(W, grid, gas, 0.0), std::invalid_argument);
    EXPECT_THROW(compute_dt(W, grid, gas, 1.5), std::invalid_argument);
}

TEST(Solver, FreeStreamPreservation) {
    const Grid grid(50, 0.0, 1.0);
    const Primitive W0{1.0, 0.7, 1.0};
    const Conserved U0 = gas.to_conserved(W0);
    for (BoundaryType bc : {BoundaryType::Periodic, BoundaryType::Transmissive}) {
        const SolverConfig cfg{FluxScheme::HLLC,
                               {ReconstructionScheme::MUSCL, Limiter::VanLeer},
                               bc, bc, 0.8, 0.5, 100000};
        Solver solver(grid, gas, cfg, std::vector<Conserved>(grid.size(), U0));
        solver.run();
        for (int j = grid.begin(); j < grid.end(); ++j) {
            EXPECT_NEAR(solver.state()[j].rho, U0.rho, 1e-13);
            EXPECT_NEAR(solver.state()[j].mom, U0.mom, 1e-13);
            EXPECT_NEAR(solver.state()[j].E, U0.E, 1e-13);
        }
    }
}

TEST(Solver, PeriodicConservation) {
    const Grid grid(64, 0.0, 1.0);
    const FluxScheme fluxes[] = {FluxScheme::Rusanov, FluxScheme::HLLC};
    const ReconstructionConfig recons[] = {
        {ReconstructionScheme::FirstOrder, Limiter::Minmod},
        {ReconstructionScheme::MUSCL, Limiter::Minmod},
        {ReconstructionScheme::MUSCL, Limiter::VanLeer}};
    for (FluxScheme fl : fluxes) {
        for (const ReconstructionConfig& rc : recons) {
            const SolverConfig cfg{fl, rc, BoundaryType::Periodic, BoundaryType::Periodic,
                                   0.8, 0.3, 100000};
            Solver solver(grid, gas, cfg, make_state(grid, sine_wave));
            const auto before = integral(solver);
            solver.run();
            const auto after = integral(solver);
            for (int q = 0; q < 3; ++q) {
                EXPECT_LE(std::abs(after[q] - before[q]), 1e-12 * std::abs(before[q]));
            }
        }
    }
}

TEST(Solver, FinalTimeAndObserver) {
    const Grid grid(64, 0.0, 1.0);
    const SolverConfig cfg{FluxScheme::HLLC,
                           {ReconstructionScheme::MUSCL, Limiter::VanLeer},
                           BoundaryType::Periodic, BoundaryType::Periodic, 0.8, 0.3, 100000};
    Solver solver(grid, gas, cfg, make_state(grid, sine_wave));
    int calls = 0;
    double last_t = 0.0;
    bool increasing = true;
    const RunStats stats = solver.run([&](const StepInfo& info) {
        ++calls;
        if (!(info.t > last_t)) increasing = false;
        last_t = info.t;
    });
    EXPECT_NEAR(stats.t_final_reached, 0.3, 1e-14);
    EXPECT_EQ(calls, stats.n_steps);
    EXPECT_TRUE(increasing);
}

TEST(Solver, SodShockTube) {
    const Grid grid(100, 0.0, 1.0);
    Solver solver(grid, gas, sod_config(), make_state(grid, sod));
    const auto before = integral(solver);
    EXPECT_NO_THROW(solver.run());
    const auto after = integral(solver);
    for (int j = grid.begin(); j < grid.end(); ++j) {
        EXPECT_TRUE(gas.is_admissible(solver.state()[j]));
    }
    EXPECT_LE(std::abs(after[0] - before[0]), 1e-12 * before[0]);
}

TEST(Solver, InvalidConfiguration) {
    const Grid grid(20, 0.0, 1.0);
    const auto U0 = make_state(grid, sod);

    SolverConfig cfg = sod_config();
    cfg.left = BoundaryType::Periodic;
    EXPECT_THROW(Solver(grid, gas, cfg, U0), std::invalid_argument);

    cfg = sod_config();
    cfg.cfl = 0.0;
    EXPECT_THROW(Solver(grid, gas, cfg, U0), std::invalid_argument);

    cfg = sod_config();
    cfg.t_final = 0.0;
    EXPECT_THROW(Solver(grid, gas, cfg, U0), std::invalid_argument);

    cfg = sod_config();
    cfg.max_steps = 0;
    EXPECT_THROW(Solver(grid, gas, cfg, U0), std::invalid_argument);

    std::vector<Conserved> wrong(grid.size() + 1);
    EXPECT_THROW(Solver(grid, gas, sod_config(), wrong), std::invalid_argument);
}

TEST(Solver, MaxStepsStop) {
    const Grid grid(100, 0.0, 1.0);
    Solver solver(grid, gas, sod_config(3), make_state(grid, sod));
    const RunStats stats = solver.run();
    EXPECT_EQ(stats.n_steps, 3);
    EXPECT_LT(stats.t_final_reached, 0.2);
}
