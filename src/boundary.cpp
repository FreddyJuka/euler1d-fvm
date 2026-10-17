#include "euler1d/boundary.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace euler1d {

BoundaryType parse_boundary_type(const std::string& name) {
    std::string lower(name);
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lower == "periodic") return BoundaryType::Periodic;
    if (lower == "transmissive") return BoundaryType::Transmissive;
    throw std::invalid_argument("Unknown boundary type '" + name +
                                "'. Valid names: periodic, transmissive");
}

void validate_boundaries(BoundaryType left, BoundaryType right) {
    if ((left == BoundaryType::Periodic) != (right == BoundaryType::Periodic)) {
        throw std::invalid_argument(
            "Periodic boundary must be applied on both sides or neither");
    }
}

void apply_boundaries(std::vector<Conserved>& U, const Grid& grid, BoundaryType left,
                      BoundaryType right) {
    if (U.size() != static_cast<std::size_t>(grid.size())) {
        throw std::invalid_argument("apply_boundaries: U.size() != grid.size()");
    }
    const int N = grid.n_cells();
    const int g = grid.n_ghost();
    for (int k = 0; k < g; ++k) {
        U[g - 1 - k] = (left == BoundaryType::Periodic) ? U[g + N - 1 - k] : U[g];
        U[g + N + k] = (right == BoundaryType::Periodic) ? U[g + k] : U[g + N - 1];
    }
}

}  // namespace euler1d
