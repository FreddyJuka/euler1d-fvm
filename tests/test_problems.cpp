#include <cmath>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "euler1d/problems.hpp"

using namespace euler1d;

namespace {

const IdealGas gas(1.4);
const double kPi = 3.14159265358979323846;

ProblemConfig wave_problem(const DensityWaveParams& w) {
    ProblemConfig p{};
    p.type = ProblemType::DensityWave;
    p.density_wave = w;
    return p;
}

}  // namespace

TEST(Problems, DensityWaveCellAverages) {
    const Grid grid(20, 0.0, 2.0);
    const DensityWaveParams w{1.0, 0.2, 1.5, 1.0, 2};
    const auto U = make_initial_condition(wave_problem(w), grid, gas);
    ASSERT_EQ(U.size(), static_cast<std::size_t>(grid.size()));
    const double kappa = 2.0 * kPi * 2 / 2.0;
    const double half = kappa * grid.dx() / 2.0;
    for (int i = grid.begin(); i < grid.end(); ++i) {
        const double rho = 1.0 + 0.2 * std::sin(kappa * grid.x_center(i)) * std::sin(half) / half;
        EXPECT_NEAR(U[i].rho, rho, 1e-14);
        EXPECT_NEAR(U[i].mom, 1.5 * rho, 1e-14);
        EXPECT_NEAR(U[i].E, 1.0 / 0.4 + 0.5 * 1.5 * 1.5 * rho, 1e-14);
    }
}

TEST(Problems, DensityWaveZeroErrorAtT0) {
    const Grid grid(32, 0.0, 1.0);
    const DensityWaveParams w{1.0, 0.2, 1.0, 1.0, 1};
    const auto U = make_initial_condition(wave_problem(w), grid, gas);
    const ErrorNorms e = density_wave_errors(U, grid, w, 0.0);
    EXPECT_NEAR(e.L1, 0.0, 1e-14);
    EXPECT_NEAR(e.L2, 0.0, 1e-14);
    EXPECT_NEAR(e.Linf, 0.0, 1e-14);
}

TEST(Problems, DensityWaveMassOverOnePeriod) {
    const Grid grid(40, -1.0, 3.0);
    const DensityWaveParams w{1.3, 0.25, 0.5, 2.0, 3};
    const auto U = make_initial_condition(wave_problem(w), grid, gas);
    double mass = 0.0;
    for (int i = grid.begin(); i < grid.end(); ++i) mass += U[i].rho * grid.dx();
    EXPECT_NEAR(mass, 1.3 * 4.0, 1e-13);
}

TEST(Problems, RiemannPlacesStates) {
    const Grid grid(10, 0.0, 1.0);
    ProblemConfig p{};
    p.type = ProblemType::Riemann;
    p.riemann = {0.5, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1};
    const auto U = make_initial_condition(p, grid, gas);
    const Conserved UL = gas.to_conserved({1.0, 0.0, 1.0});
    const Conserved UR = gas.to_conserved({0.125, 0.0, 0.1});
    for (int i = grid.begin(); i < grid.end(); ++i) {
        const Conserved& ref = grid.x_center(i) < 0.5 ? UL : UR;
        EXPECT_EQ(U[i].rho, ref.rho);
        EXPECT_EQ(U[i].mom, ref.mom);
        EXPECT_EQ(U[i].E, ref.E);
    }
}

TEST(Problems, InvalidDensityWaveThrows) {
    const Grid grid(10, 0.0, 1.0);
    EXPECT_THROW(make_initial_condition(wave_problem({1.0, 1.0, 1.0, 1.0, 1}), grid, gas),
                 std::invalid_argument);
    EXPECT_THROW(make_initial_condition(wave_problem({1.0, 0.2, 1.0, 0.0, 1}), grid, gas),
                 std::invalid_argument);
    EXPECT_THROW(make_initial_condition(wave_problem({1.0, 0.2, 1.0, 1.0, 0}), grid, gas),
                 std::invalid_argument);
}

TEST(Problems, ParseProblemType) {
    EXPECT_EQ(parse_problem_type("Riemann"), ProblemType::Riemann);
    EXPECT_EQ(parse_problem_type("DENSITY_WAVE"), ProblemType::DensityWave);
    EXPECT_THROW(parse_problem_type("sedov"), std::invalid_argument);
}

TEST(Problems, InvalidRiemannThrows) {
    ProblemConfig p{};
    p.type = ProblemType::Riemann;
    p.riemann = {0.5, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1};
    const Grid grid(10, 0.0, 1.0);
    EXPECT_NO_THROW(make_initial_condition(p, grid, gas));

    ProblemConfig bad = p;
    bad.riemann.rho_L = -1.0;
    EXPECT_THROW(make_initial_condition(bad, grid, gas), std::invalid_argument);
    bad = p;
    bad.riemann.p_R = 0.0;
    EXPECT_THROW(make_initial_condition(bad, grid, gas), std::invalid_argument);
    bad = p;
    bad.riemann.u_L = std::nan("");
    EXPECT_THROW(make_initial_condition(bad, grid, gas), std::invalid_argument);
}
