#include <algorithm>
#include <random>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "euler1d/mesh.hpp"
#include "euler1d/reconstruction.hpp"
#include "euler1d/types.hpp"

using namespace euler1d;

namespace {

const Limiter kLimiters[] = {Limiter::Minmod, Limiter::VanLeer};
const Limiter kLinearLimiters[] = {Limiter::Minmod, Limiter::VanLeer, Limiter::Unlimited};

void expect_prim_eq(const Primitive& a, const Primitive& b) {
    EXPECT_EQ(a.rho, b.rho);
    EXPECT_EQ(a.u, b.u);
    EXPECT_EQ(a.p, b.p);
}

void expect_prim_near(const Primitive& a, const Primitive& b, double tol) {
    EXPECT_NEAR(a.rho, b.rho, tol);
    EXPECT_NEAR(a.u, b.u, tol);
    EXPECT_NEAR(a.p, b.p, tol);
}

void expect_within(double v, double a, double b) {
    EXPECT_GE(v, std::min(a, b) - 1e-14);
    EXPECT_LE(v, std::max(a, b) + 1e-14);
}

}  // namespace

TEST(Limiters, Values) {
    EXPECT_DOUBLE_EQ(minmod(1.0, 2.0), 1.0);
    EXPECT_DOUBLE_EQ(minmod(-1.0, -3.0), -1.0);
    EXPECT_DOUBLE_EQ(minmod(1.0, -1.0), 0.0);
    EXPECT_DOUBLE_EQ(minmod(0.0, 1.0), 0.0);
    EXPECT_DOUBLE_EQ(van_leer(1.0, 1.0), 1.0);
    EXPECT_DOUBLE_EQ(van_leer(1.0, 3.0), 1.5);
    EXPECT_DOUBLE_EQ(van_leer(1.0, -1.0), 0.0);
    EXPECT_DOUBLE_EQ(van_leer(0.0, 2.0), 0.0);
    EXPECT_DOUBLE_EQ(unlimited_slope(1.0, 3.0), 2.0);
    EXPECT_DOUBLE_EQ(unlimited_slope(1.0, -1.0), 0.0);
    EXPECT_DOUBLE_EQ(unlimited_slope(-2.0, -4.0), -3.0);
}

TEST(Reconstruction, FirstOrder) {
    const Grid grid(6, 0.0, 1.0);
    std::vector<Primitive> W(grid.size());
    for (int i = 0; i < grid.size(); ++i) W[i] = {1.0 + i, 0.1 * i, 2.0 + 3.0 * i};
    std::vector<Primitive> WL, WR;
    reconstruct(W, grid, {ReconstructionScheme::FirstOrder, Limiter::Minmod}, WL, WR);
    ASSERT_EQ(WL.size(), 7u);
    ASSERT_EQ(WR.size(), 7u);
    for (int f = 0; f <= 6; ++f) {
        expect_prim_eq(WL[f], W[grid.begin() + f - 1]);
        expect_prim_eq(WR[f], W[grid.begin() + f]);
    }
}

TEST(Reconstruction, LinearExactness) {
    const Grid grid(8, 0.0, 1.0);
    std::vector<Primitive> W(grid.size());
    for (int i = 0; i < grid.size(); ++i) W[i] = {1.0 + 0.1 * i, 0.5 + 0.01 * i, 2.0 + 0.2 * i};
    for (Limiter lim : kLinearLimiters) {
        std::vector<Primitive> WL, WR;
        reconstruct(W, grid, {ReconstructionScheme::MUSCL, lim}, WL, WR);
        for (int f = 0; f <= grid.n_cells(); ++f) {
            const int i = grid.begin() + f;
            const Primitive mid{0.5 * (W[i - 1].rho + W[i].rho), 0.5 * (W[i - 1].u + W[i].u),
                                0.5 * (W[i - 1].p + W[i].p)};
            expect_prim_near(WL[f], mid, 1e-13);
            expect_prim_near(WR[f], mid, 1e-13);
        }
    }
}

TEST(Reconstruction, StepPreservation) {
    const Grid grid(8, 0.0, 1.0);
    std::vector<Primitive> W(grid.size());
    for (int i = 0; i < grid.size(); ++i) {
        W[i] = (i < 6) ? Primitive{1.0, 0.0, 1.0} : Primitive{0.125, 0.0, 0.1};
    }
    for (Limiter lim : kLimiters) {
        std::vector<Primitive> WL, WR;
        reconstruct(W, grid, {ReconstructionScheme::MUSCL, lim}, WL, WR);
        const int f = 6 - grid.begin();
        expect_prim_eq(WL[f], W[5]);
        expect_prim_eq(WR[f], W[6]);
    }
}

TEST(Reconstruction, UnlimitedIsNotBounded) {
    const Grid grid(8, 0.0, 1.0);
    std::vector<Primitive> W(grid.size());
    for (int i = 0; i < grid.size(); ++i) {
        W[i] = (i < 6) ? Primitive{1.0, 0.0, 1.0} : Primitive{0.125, 0.0, 0.1};
    }
    std::vector<Primitive> WL, WR;
    reconstruct(W, grid, {ReconstructionScheme::MUSCL, Limiter::Unlimited}, WL, WR);
    bool outside = false;
    for (int f = 0; f <= grid.n_cells(); ++f) {
        const Primitive& a = W[grid.begin() + f - 1];
        const Primitive& b = W[grid.begin() + f];
        for (const Primitive* v : {&WL[f], &WR[f]}) {
            if (v->rho < std::min(a.rho, b.rho) - 1e-14 ||
                v->rho > std::max(a.rho, b.rho) + 1e-14) {
                outside = true;
            }
        }
    }
    EXPECT_TRUE(outside);
}

TEST(Reconstruction, Boundedness) {
    const Grid grid(64, 0.0, 1.0);
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> pos(0.1, 10.0), vel(-2.0, 2.0);
    std::vector<Primitive> W(grid.size());
    for (auto& w : W) w = {pos(rng), vel(rng), pos(rng)};
    for (Limiter lim : kLimiters) {
        std::vector<Primitive> WL, WR;
        reconstruct(W, grid, {ReconstructionScheme::MUSCL, lim}, WL, WR);
        for (int f = 0; f <= grid.n_cells(); ++f) {
            const Primitive& a = W[grid.begin() + f - 1];
            const Primitive& b = W[grid.begin() + f];
            for (const Primitive* v : {&WL[f], &WR[f]}) {
                expect_within(v->rho, a.rho, b.rho);
                expect_within(v->u, a.u, b.u);
                expect_within(v->p, a.p, b.p);
            }
        }
    }
}

TEST(Reconstruction, InvalidInputs) {
    const ReconstructionConfig muscl_cfg{ReconstructionScheme::MUSCL, Limiter::Minmod};
    const Grid thin(6, 0.0, 1.0, 1);
    std::vector<Primitive> W(thin.size(), Primitive{1.0, 0.0, 1.0});
    std::vector<Primitive> WL, WR;
    EXPECT_THROW(reconstruct(W, thin, muscl_cfg, WL, WR), std::invalid_argument);

    const Grid grid(6, 0.0, 1.0);
    std::vector<Primitive> bad(grid.size() + 1, Primitive{1.0, 0.0, 1.0});
    EXPECT_THROW(reconstruct(bad, grid, muscl_cfg, WL, WR), std::invalid_argument);
}

TEST(Reconstruction, Parsers) {
    EXPECT_EQ(parse_reconstruction("MUSCL"), ReconstructionScheme::MUSCL);
    EXPECT_EQ(parse_reconstruction("First_Order"), ReconstructionScheme::FirstOrder);
    EXPECT_EQ(parse_limiter("MinMod"), Limiter::Minmod);
    EXPECT_EQ(parse_limiter("VanLeer"), Limiter::VanLeer);
    EXPECT_EQ(parse_limiter("van_leer"), Limiter::VanLeer);
    EXPECT_EQ(parse_limiter("Unlimited"), Limiter::Unlimited);
    EXPECT_THROW(parse_reconstruction("weno"), std::invalid_argument);
    EXPECT_THROW(parse_limiter("superbee"), std::invalid_argument);
}
