#pragma once

// Input and output: reading case configuration files and writing solution
// data for post-processing.

#include <fstream>
#include <istream>
#include <optional>
#include <string>
#include <vector>

#include "euler1d/eos.hpp"
#include "euler1d/mesh.hpp"
#include "euler1d/problems.hpp"
#include "euler1d/solver.hpp"
#include "euler1d/types.hpp"

namespace euler1d {

struct Config {
    ProblemConfig problem;
    int n_cells;
    double x_min;
    double x_max;
    double gamma;
    SolverConfig solver;
    std::string problem_name;
    std::string flux_name;
    std::string reconstruction_name;
    std::string limiter_name;  // empty if not given
    std::string output_dir;
};

// Parses "key = value" lines ('#' starts a comment). Throws std::runtime_error with the
// line number on malformed lines, duplicate/unknown keys, missing keys or bad values.
Config parse_config(std::istream& in);
Config load_config(const std::string& path);

void write_solution_csv(const std::string& path, const std::vector<Conserved>& U,
                        const Grid& grid, const IdealGas& gas);

class ResidualWriter {
public:
    explicit ResidualWriter(const std::string& path);
    void write(const StepInfo& info);

private:
    std::ofstream out_;
};

struct SummaryInfo {
    std::string problem;
    int n_cells;
    double dx;
    std::string flux;
    std::string reconstruction;
    std::string limiter;  // "none" if unused
    double cfl;
    double t_final;
    double t_reached;
    int n_steps;
    double wall_time_s;
    std::optional<ErrorNorms> rho_errors;  // density_wave only
};

void write_summary_csv(const std::string& path, const SummaryInfo& info);

}  // namespace euler1d
