#include <cmath>
#include <stdexcept>

#include <gtest/gtest.h>

#include "euler1d/eos.hpp"
#include "euler1d/types.hpp"

using namespace euler1d;

namespace {

void expect_close(double actual, double exact) {
    if (exact == 0.0) {
        EXPECT_LT(std::abs(actual), 1e-14);
    } else {
        EXPECT_LT(std::abs(actual - exact) / std::abs(exact), 1e-13);
    }
}

}  // namespace

TEST(IdealGas, ConstructorRejectsGammaAtMostOne) {
    EXPECT_THROW(IdealGas(1.0), std::invalid_argument);
    EXPECT_THROW(IdealGas(0.5), std::invalid_argument);
    EXPECT_NO_THROW(IdealGas(1.4));
}

TEST(IdealGas, ToConservedAtRest) {
    const IdealGas eos(1.4);
    const Conserved U = eos.to_conserved({1.0, 0.0, 1.0});
    EXPECT_DOUBLE_EQ(U.rho, 1.0);
    EXPECT_DOUBLE_EQ(U.mom, 0.0);
    EXPECT_DOUBLE_EQ(U.E, 2.5);
}

TEST(IdealGas, ToConservedMoving) {
    const IdealGas eos(1.4);
    const Conserved U = eos.to_conserved({1.0, 2.0, 1.0});
    EXPECT_DOUBLE_EQ(U.rho, 1.0);
    EXPECT_DOUBLE_EQ(U.mom, 2.0);
    EXPECT_DOUBLE_EQ(U.E, 4.5);
}

TEST(IdealGas, RoundTrip) {
    const IdealGas eos(1.4);
    const Primitive states[] = {
        {0.125, -0.75, 0.1}, {1.0, 2.0, 1.0}, {5.99924, 19.5975, 460.894}};
    for (const Primitive& W : states) {
        const Primitive R = eos.to_primitive(eos.to_conserved(W));
        expect_close(R.rho, W.rho);
        expect_close(R.u, W.u);
        expect_close(R.p, W.p);
    }
}

TEST(IdealGas, SoundSpeed) {
    const IdealGas eos(1.4);
    EXPECT_NEAR(eos.sound_speed({1.0, 0.0, 1.0}), std::sqrt(1.4), 1e-14);
}

TEST(IdealGas, Admissibility) {
    const IdealGas eos(1.4);
    EXPECT_TRUE(eos.is_admissible({1.0, 0.0, 2.5}));
    EXPECT_FALSE(eos.is_admissible({1.0, 2.0, 1.0}));
    EXPECT_FALSE(eos.is_admissible({-1.0, 0.0, 2.5}));
}

TEST(Conserved, Arithmetic) {
    const Conserved U1{1.5, -2.0, 4.25};
    const Conserved U2{0.5, 3.0, 1.75};

    const Conserved s = U1 + U2;
    EXPECT_EQ(s.rho, 2.0);
    EXPECT_EQ(s.mom, 1.0);
    EXPECT_EQ(s.E, 6.0);

    const Conserved d = U1 - U2;
    EXPECT_EQ(d.rho, 1.0);
    EXPECT_EQ(d.mom, -5.0);
    EXPECT_EQ(d.E, 2.5);

    for (const Conserved& m : {2.0 * U1, U1 * 2.0}) {
        EXPECT_EQ(m.rho, 3.0);
        EXPECT_EQ(m.mom, -4.0);
        EXPECT_EQ(m.E, 8.5);
    }

    Conserved a = U1;
    a += U2;
    EXPECT_EQ(a.E, 6.0);
    a -= U2;
    EXPECT_EQ(a.E, 4.25);
}
