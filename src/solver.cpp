#include "euler1d/solver.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "euler1d/time_integration.hpp"

namespace euler1d {

Solver::Solver(const Grid& grid, const IdealGas& gas, const SolverConfig& cfg,
               std::vector<Conserved> U0)
    : grid_(grid), gas_(gas), cfg_(cfg), flux_(select_flux(cfg.flux)), U_(std::move(U0)) {
    validate_boundaries(cfg.left, cfg.right);
    if (!(cfg.cfl > 0.0 && cfg.cfl <= 1.0)) {
        throw std::invalid_argument("Solver: cfl must be in (0, 1]");
    }
    if (!(cfg.t_final > 0.0)) throw std::invalid_argument("Solver: t_final must be > 0");
    if (cfg.max_steps <= 0) throw std::invalid_argument("Solver: max_steps must be > 0");
    if (U_.size() != static_cast<std::size_t>(grid.size())) {
        throw std::invalid_argument("Solver: U0.size() != grid.size()");
    }
    if (cfg.recon.scheme == ReconstructionScheme::MUSCL && grid.n_ghost() < 2) {
        throw std::invalid_argument("Solver: MUSCL requires n_ghost >= 2");
    }

    const std::size_t n = U_.size();
    const std::size_t n_faces = static_cast<std::size_t>(grid.n_cells()) + 1;
    U1_.resize(n);
    U2_.resize(n);
    L_.resize(n);
    Uw_.resize(n);
    W_.resize(n);
    WL_.resize(n_faces);
    WR_.resize(n_faces);
    F_.resize(n_faces);
}

void Solver::spatial_operator(const std::vector<Conserved>& U_in, std::vector<Conserved>& L) {
    Uw_ = U_in;
    apply_boundaries(Uw_, grid_, cfg_.left, cfg_.right);
    for (int i = 0; i < grid_.size(); ++i) W_[i] = gas_.to_primitive(Uw_[i]);
    reconstruct(W_, grid_, cfg_.recon, WL_, WR_);

    const int N = grid_.n_cells();
    for (int f = 0; f <= N; ++f) F_[f] = flux_(WL_[f], WR_[f], gas_);

    const double inv_dx = 1.0 / grid_.dx();
    for (int j = grid_.begin(); j < grid_.end(); ++j) {
        const int f = j - grid_.begin();
        L[j] = (-inv_dx) * (F_[f + 1] - F_[f]);
    }
}

void Solver::check_admissible() const {
    for (int j = grid_.begin(); j < grid_.end(); ++j) {
        if (!gas_.is_admissible(U_[j])) {
            const Primitive W = gas_.to_primitive(U_[j]);
            std::ostringstream msg;
            msg << "Solver: inadmissible state at step " << step_count_ << ", t = " << t_
                << ", cell index " << j << ": rho = " << W.rho << ", p = " << W.p;
            throw std::runtime_error(msg.str());
        }
    }
}

double Solver::step(double dt_max) {
    for (int j = grid_.begin(); j < grid_.end(); ++j) W_[j] = gas_.to_primitive(U_[j]);
    const double dt = std::min(compute_dt(W_, grid_, gas_, cfg_.cfl), dt_max);

    const int b = grid_.begin();
    const int e = grid_.end();

    // Stage 1: U1 = Un + dt*L(Un)
    spatial_operator(U_, L_);
    for (int j = b; j < e; ++j) U1_[j] = U_[j] + dt * L_[j];

    // Stage 2: U2 = 3/4 Un + 1/4 (U1 + dt*L(U1))
    spatial_operator(U1_, L_);
    for (int j = b; j < e; ++j) U2_[j] = 0.75 * U_[j] + 0.25 * (U1_[j] + dt * L_[j]);

    // Stage 3: Un+1 = 1/3 Un + 2/3 (U2 + dt*L(U2)); residual uses Un before overwrite.
    spatial_operator(U2_, L_);
    std::array<double, 3> sum{0.0, 0.0, 0.0};
    for (int j = b; j < e; ++j) {
        const Conserved Unew = (1.0 / 3.0) * U_[j] + (2.0 / 3.0) * (U2_[j] + dt * L_[j]);
        const Conserved r = (1.0 / dt) * (Unew - U_[j]);
        sum[0] += r.rho * r.rho;
        sum[1] += r.mom * r.mom;
        sum[2] += r.E * r.E;
        U_[j] = Unew;
    }
    for (int q = 0; q < 3; ++q) last_residual_[q] = std::sqrt(sum[q] * grid_.dx());

    t_ += dt;
    ++step_count_;
    check_admissible();
    return dt;
}

RunStats Solver::run(const std::function<void(const StepInfo&)>& observer) {
    int n_steps = 0;
    while (t_ < cfg_.t_final && n_steps < cfg_.max_steps) {
        const double dt_max = cfg_.t_final - t_;
        const double dt = step(dt_max);
        if (dt == dt_max) t_ = cfg_.t_final;
        ++n_steps;
        if (observer) observer({step_count_, t_, dt, last_residual_});
    }
    return {n_steps, t_};
}

}  // namespace euler1d
