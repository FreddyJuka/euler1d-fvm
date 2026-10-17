#pragma once

// Piecewise-constant and MUSCL reconstruction of face states with slope limiters.

#include <string>
#include <vector>

#include "euler1d/mesh.hpp"
#include "euler1d/types.hpp"

namespace euler1d {

enum class ReconstructionScheme { FirstOrder, MUSCL };
enum class Limiter { Minmod, VanLeer };

struct ReconstructionConfig {
    ReconstructionScheme scheme;
    Limiter limiter;
};

// 0 if a*b <= 0, else sign(a)*min(|a|,|b|).
double minmod(double a, double b);

// 0 if a*b <= 0, else 2*a*b/(a+b).
double van_leer(double a, double b);

// Case-insensitive; throw std::invalid_argument for unknown names.
ReconstructionScheme parse_reconstruction(const std::string& name);
Limiter parse_limiter(const std::string& name);

// Reconstructs left/right states at the N+1 faces; face f lies between storage
// cells i-1 and i with i = grid.begin() + f. W must include filled ghost cells.
// Throws std::invalid_argument if W.size() != grid.size(), or for MUSCL if
// grid.n_ghost() < 2. WL and WR are resized only if their size differs.
void reconstruct(const std::vector<Primitive>& W, const Grid& grid,
                 const ReconstructionConfig& cfg, std::vector<Primitive>& WL,
                 std::vector<Primitive>& WR);

}  // namespace euler1d
