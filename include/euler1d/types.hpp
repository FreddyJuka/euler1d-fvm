#pragma once

// Core data structures shared across the solver: conservative and primitive
// state vectors.

namespace euler1d {

// Conservative variables: density, momentum (rho*u), total energy per volume.
struct Conserved {
    double rho;
    double mom;
    double E;
};

// Primitive variables: density, velocity, pressure.
struct Primitive {
    double rho;
    double u;
    double p;
};

constexpr Conserved operator+(const Conserved& a, const Conserved& b) {
    return {a.rho + b.rho, a.mom + b.mom, a.E + b.E};
}

constexpr Conserved operator-(const Conserved& a, const Conserved& b) {
    return {a.rho - b.rho, a.mom - b.mom, a.E - b.E};
}

constexpr Conserved& operator+=(Conserved& a, const Conserved& b) {
    a.rho += b.rho;
    a.mom += b.mom;
    a.E += b.E;
    return a;
}

constexpr Conserved& operator-=(Conserved& a, const Conserved& b) {
    a.rho -= b.rho;
    a.mom -= b.mom;
    a.E -= b.E;
    return a;
}

constexpr Conserved operator*(double s, const Conserved& a) {
    return {s * a.rho, s * a.mom, s * a.E};
}

constexpr Conserved operator*(const Conserved& a, double s) {
    return s * a;
}

}  // namespace euler1d
