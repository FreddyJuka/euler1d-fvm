#pragma once

// Boundary condition application via ghost cells: periodic and transmissive.

#include <string>
#include <vector>

#include "euler1d/mesh.hpp"
#include "euler1d/types.hpp"

namespace euler1d {

enum class BoundaryType { Periodic, Transmissive };

// Case-insensitive; throws std::invalid_argument for unknown names.
BoundaryType parse_boundary_type(const std::string& name);

// Throws std::invalid_argument if exactly one side is Periodic.
void validate_boundaries(BoundaryType left, BoundaryType right);

// Fills ghost cells; interior cells are never modified.
// Throws std::invalid_argument if U.size() != grid.size().
void apply_boundaries(std::vector<Conserved>& U, const Grid& grid, BoundaryType left,
                      BoundaryType right);

}  // namespace euler1d
