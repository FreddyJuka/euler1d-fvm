#include "euler1d/reconstruction.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace euler1d {

double minmod(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return std::copysign(std::min(std::abs(a), std::abs(b)), a);
}

double van_leer(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return 2.0 * a * b / (a + b);
}

double unlimited_slope(double a, double b) { return 0.5 * (a + b); }

namespace {

std::string to_lower(const std::string& s) {
    std::string out(s);
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return out;
}

// Limited slope sigma of cell c for each primitive component.
template <typename LimiterFn>
Primitive limited_slope(const std::vector<Primitive>& W, int c, LimiterFn limit) {
    const Primitive& m = W[c - 1];
    const Primitive& o = W[c];
    const Primitive& p = W[c + 1];
    return {limit(o.rho - m.rho, p.rho - o.rho), limit(o.u - m.u, p.u - o.u),
            limit(o.p - m.p, p.p - o.p)};
}

template <typename LimiterFn>
void muscl(const std::vector<Primitive>& W, const Grid& grid, std::vector<Primitive>& WL,
           std::vector<Primitive>& WR, LimiterFn limit) {
    for (int f = 0; f <= grid.n_cells(); ++f) {
        const int i = grid.begin() + f;
        const Primitive sl = limited_slope(W, i - 1, limit);
        const Primitive sr = limited_slope(W, i, limit);
        const Primitive& a = W[i - 1];
        const Primitive& b = W[i];
        WL[f] = {a.rho + 0.5 * sl.rho, a.u + 0.5 * sl.u, a.p + 0.5 * sl.p};
        WR[f] = {b.rho - 0.5 * sr.rho, b.u - 0.5 * sr.u, b.p - 0.5 * sr.p};
    }
}

}  // namespace

ReconstructionScheme parse_reconstruction(const std::string& name) {
    const std::string lower = to_lower(name);
    if (lower == "first_order") return ReconstructionScheme::FirstOrder;
    if (lower == "muscl") return ReconstructionScheme::MUSCL;
    throw std::invalid_argument("Unknown reconstruction scheme '" + name +
                                "'. Valid names: first_order, muscl");
}

Limiter parse_limiter(const std::string& name) {
    const std::string lower = to_lower(name);
    if (lower == "minmod") return Limiter::Minmod;
    if (lower == "vanleer" || lower == "van_leer") return Limiter::VanLeer;
    if (lower == "unlimited") return Limiter::Unlimited;
    throw std::invalid_argument("Unknown limiter '" + name +
                                "'. Valid names: minmod, vanleer, van_leer, unlimited");
}

void reconstruct(const std::vector<Primitive>& W, const Grid& grid,
                 const ReconstructionConfig& cfg, std::vector<Primitive>& WL,
                 std::vector<Primitive>& WR) {
    if (W.size() != static_cast<std::size_t>(grid.size())) {
        throw std::invalid_argument("reconstruct: W.size() != grid.size()");
    }
    if (cfg.scheme == ReconstructionScheme::MUSCL && grid.n_ghost() < 2) {
        throw std::invalid_argument("reconstruct: MUSCL requires n_ghost >= 2");
    }

    const std::size_t n_faces = static_cast<std::size_t>(grid.n_cells()) + 1;
    if (WL.size() != n_faces) WL.resize(n_faces);
    if (WR.size() != n_faces) WR.resize(n_faces);

    if (cfg.scheme == ReconstructionScheme::FirstOrder) {
        for (int f = 0; f <= grid.n_cells(); ++f) {
            const int i = grid.begin() + f;
            WL[f] = W[i - 1];
            WR[f] = W[i];
        }
        return;
    }

    switch (cfg.limiter) {
        case Limiter::Minmod:
            muscl(W, grid, WL, WR, &minmod);
            break;
        case Limiter::VanLeer:
            muscl(W, grid, WL, WR, &van_leer);
            break;
        case Limiter::Unlimited:
            muscl(W, grid, WL, WR, &unlimited_slope);
            break;
    }
}

}  // namespace euler1d
