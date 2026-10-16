#include <cmath>
#include <stdexcept>

#include <gtest/gtest.h>

#include "euler1d/eos.hpp"
#include "euler1d/flux.hpp"
#include "euler1d/types.hpp"

using namespace euler1d;

namespace {

const IdealGas gas(1.4);

void expect_rel(double actual, double exact) {
    const double scale = std::max(std::abs(exact), 1.0);
    EXPECT_LE(std::abs(actual - exact), 1e-14 * scale);
}

Primitive mirror(const Primitive& W) { return {W.rho, -W.u, W.p}; }

}  // namespace

TEST(PhysicalFlux, AtRest) {
    const Conserved F = physical_flux({1.0, 0.0, 1.0}, gas);
    EXPECT_NEAR(F.rho, 0.0, 1e-14);
    EXPECT_NEAR(F.mom, 1.0, 1e-14);
    EXPECT_NEAR(F.E, 0.0, 1e-14);
}

TEST(PhysicalFlux, Moving) {
    const Conserved F = physical_flux({1.0, 2.0, 1.0}, gas);
    EXPECT_NEAR(F.rho, 2.0, 1e-14);
    EXPECT_NEAR(F.mom, 5.0, 1e-14);
    EXPECT_NEAR(F.E, 11.0, 1e-14);
}

TEST(RusanovFlux, Consistency) {
    const Primitive states[] = {
        {1.0, 0.0, 1.0}, {0.125, -0.75, 0.1}, {5.99924, 19.5975, 460.894}};
    for (const Primitive& W : states) {
        const Conserved F = rusanov_flux(W, W, gas);
        const Conserved P = physical_flux(W, gas);
        expect_rel(F.rho, P.rho);
        expect_rel(F.mom, P.mom);
        expect_rel(F.E, P.E);
    }
}

TEST(RusanovFlux, SodInterface) {
    const Conserved F = rusanov_flux({1.0, 0.0, 1.0}, {0.125, 0.0, 0.1}, gas);
    EXPECT_NEAR(F.rho, 0.5 * std::sqrt(1.4) * 0.875, 1e-13);
    EXPECT_NEAR(F.mom, 0.55, 1e-13);
    EXPECT_NEAR(F.E, 0.5 * std::sqrt(1.4) * 2.25, 1e-13);
}

TEST(RusanovFlux, MirrorSymmetry) {
    const Primitive WL{1.0, 0.3, 1.0};
    const Primitive WR{0.5, -0.2, 0.4};
    const Conserved F = rusanov_flux(WL, WR, gas);
    const Conserved G = rusanov_flux(mirror(WR), mirror(WL), gas);
    EXPECT_NEAR(F.rho, -G.rho, 1e-14);
    EXPECT_NEAR(F.mom, G.mom, 1e-14);
    EXPECT_NEAR(F.E, -G.E, 1e-14);
}

TEST(FluxScheme, Parse) {
    EXPECT_EQ(parse_flux_scheme("Rusanov"), FluxScheme::Rusanov);
    EXPECT_EQ(parse_flux_scheme("RUSANOV"), FluxScheme::Rusanov);
    EXPECT_THROW(parse_flux_scheme("roe"), std::invalid_argument);
}

TEST(FluxScheme, Select) {
    const FluxFunction f = select_flux(FluxScheme::Rusanov);
    const Primitive WL{1.0, 0.3, 1.0};
    const Primitive WR{0.5, -0.2, 0.4};
    const Conserved a = f(WL, WR, gas);
    const Conserved b = rusanov_flux(WL, WR, gas);
    EXPECT_EQ(a.rho, b.rho);
    EXPECT_EQ(a.mom, b.mom);
    EXPECT_EQ(a.E, b.E);
}
