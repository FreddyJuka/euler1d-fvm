#include <cmath>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "euler1d/boundary.hpp"
#include "euler1d/mesh.hpp"
#include "euler1d/types.hpp"

using namespace euler1d;

namespace {

Conserved fill_value(int j) { return {j + 1.0, 10.0 * (j + 1), 100.0 * (j + 1)}; }

std::vector<Conserved> make_filled(const Grid& grid) {
    std::vector<Conserved> U(grid.size(), Conserved{-1.0, -1.0, -1.0});
    for (int j = 0; j < grid.n_cells(); ++j) U[grid.begin() + j] = fill_value(j);
    return U;
}

void expect_equal(const Conserved& a, const Conserved& b) {
    EXPECT_EQ(a.rho, b.rho);
    EXPECT_EQ(a.mom, b.mom);
    EXPECT_EQ(a.E, b.E);
}

void expect_interior_unchanged(const std::vector<Conserved>& U, const Grid& grid) {
    for (int j = 0; j < grid.n_cells(); ++j) expect_equal(U[grid.begin() + j], fill_value(j));
}

}  // namespace

TEST(Grid, Geometry) {
    const Grid grid(10, 0.0, 1.0);
    EXPECT_NEAR(grid.dx(), 0.1, 1e-14);
    EXPECT_EQ(grid.size(), 14);
    EXPECT_EQ(grid.begin(), 2);
    EXPECT_EQ(grid.end(), 12);
    EXPECT_NEAR(grid.x_center(2), 0.05, 1e-14);
    EXPECT_NEAR(grid.x_center(11), 0.95, 1e-14);
    EXPECT_NEAR(grid.x_center(0), -0.15, 1e-14);
    EXPECT_NEAR(grid.x_center(13), 1.15, 1e-14);
}

TEST(Grid, InvalidArguments) {
    EXPECT_THROW(Grid(0, 0.0, 1.0), std::invalid_argument);
    EXPECT_THROW(Grid(10, 1.0, 1.0), std::invalid_argument);
    EXPECT_THROW(Grid(10, 0.0, 1.0, 0), std::invalid_argument);
}

TEST(Boundary, Periodic) {
    const Grid grid(5, 0.0, 1.0);
    auto U = make_filled(grid);
    apply_boundaries(U, grid, BoundaryType::Periodic, BoundaryType::Periodic);
    expect_equal(U[0], fill_value(3));
    expect_equal(U[1], fill_value(4));
    expect_equal(U[7], fill_value(0));
    expect_equal(U[8], fill_value(1));
    expect_interior_unchanged(U, grid);
}

TEST(Boundary, Transmissive) {
    const Grid grid(5, 0.0, 1.0);
    auto U = make_filled(grid);
    apply_boundaries(U, grid, BoundaryType::Transmissive, BoundaryType::Transmissive);
    expect_equal(U[0], fill_value(0));
    expect_equal(U[1], fill_value(0));
    expect_equal(U[7], fill_value(4));
    expect_equal(U[8], fill_value(4));
    expect_interior_unchanged(U, grid);
}

TEST(Boundary, TransmissiveThreeGhosts) {
    const Grid grid(5, 0.0, 1.0, 3);
    auto U = make_filled(grid);
    apply_boundaries(U, grid, BoundaryType::Transmissive, BoundaryType::Transmissive);
    for (int k = 0; k < 3; ++k) {
        expect_equal(U[k], fill_value(0));
        expect_equal(U[grid.end() + k], fill_value(4));
    }
    expect_interior_unchanged(U, grid);
}

TEST(Boundary, ValidateBoundaries) {
    EXPECT_THROW(validate_boundaries(BoundaryType::Periodic, BoundaryType::Transmissive),
                 std::invalid_argument);
    EXPECT_THROW(validate_boundaries(BoundaryType::Transmissive, BoundaryType::Periodic),
                 std::invalid_argument);
    EXPECT_NO_THROW(validate_boundaries(BoundaryType::Periodic, BoundaryType::Periodic));
    EXPECT_NO_THROW(
        validate_boundaries(BoundaryType::Transmissive, BoundaryType::Transmissive));
}

TEST(Boundary, SizeMismatchThrows) {
    const Grid grid(5, 0.0, 1.0);
    std::vector<Conserved> U(grid.size() + 1);
    EXPECT_THROW(
        apply_boundaries(U, grid, BoundaryType::Transmissive, BoundaryType::Transmissive),
        std::invalid_argument);
}

TEST(Boundary, Parse) {
    EXPECT_EQ(parse_boundary_type("Periodic"), BoundaryType::Periodic);
    EXPECT_EQ(parse_boundary_type("TRANSMISSIVE"), BoundaryType::Transmissive);
    EXPECT_THROW(parse_boundary_type("reflective"), std::invalid_argument);
}
