#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

#include "euler1d/io.hpp"
#include "euler1d/solver.hpp"

using namespace euler1d;

namespace {

void print_usage(std::ostream& os) {
    os << "Usage: euler1d <config.cfg> [--n-cells N] [--output-dir DIR] [--flux NAME]\n"
          "                [--reconstruction NAME] [--limiter NAME] [--quiet] [--help]\n"
          "  --n-cells N            override n_cells from the config file\n"
          "  --output-dir DIR       override output_dir from the config file\n"
          "  --flux NAME            override flux (rusanov, hllc)\n"
          "  --reconstruction NAME  override reconstruction (first_order, muscl)\n"
          "  --limiter NAME         override limiter (minmod, vanleer, van_leer)\n"
          "  --quiet                suppress the final summary line\n"
          "  --help                 show this message\n";
}

int parse_int_arg(const std::string& text) {
    std::size_t pos = 0;
    int v = 0;
    try {
        v = std::stoi(text, &pos);
    } catch (const std::exception&) {
        pos = 0;
    }
    if (pos == 0 || pos != text.size()) {
        throw std::runtime_error("invalid integer '" + text + "' for --n-cells");
    }
    return v;
}

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

int run_cli(int argc, char** argv) {
    std::string config_path;
    std::optional<int> n_cells;
    std::optional<std::string> output_dir;
    std::optional<std::string> flux, reconstruction, limiter;
    bool quiet = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help") {
            print_usage(std::cout);
            return 0;
        } else if (arg == "--quiet") {
            quiet = true;
        } else if (arg == "--n-cells" || arg == "--output-dir" || arg == "--flux" ||
                   arg == "--reconstruction" || arg == "--limiter") {
            if (i + 1 >= argc) throw std::runtime_error("option " + arg + " requires a value");
            const std::string value = argv[++i];
            if (arg == "--n-cells") n_cells = parse_int_arg(value);
            else if (arg == "--output-dir") output_dir = value;
            else if (arg == "--flux") flux = value;
            else if (arg == "--reconstruction") reconstruction = value;
            else limiter = value;
        } else if (!arg.empty() && arg[0] == '-') {
            throw std::runtime_error("unknown option '" + arg + "'");
        } else if (config_path.empty()) {
            config_path = arg;
        } else {
            throw std::runtime_error("unexpected argument '" + arg + "'");
        }
    }
    if (config_path.empty()) throw std::runtime_error("missing config file (see --help)");

    Config cfg = load_config(config_path);
    if (n_cells) cfg.n_cells = *n_cells;
    if (output_dir) cfg.output_dir = *output_dir;
    if (flux) {
        cfg.solver.flux = parse_flux_scheme(*flux);
        cfg.flux_name = to_lower(*flux);
    }
    if (reconstruction) {
        cfg.solver.recon.scheme = parse_reconstruction(*reconstruction);
        cfg.reconstruction_name = to_lower(*reconstruction);
    }
    if (limiter) {
        cfg.solver.recon.limiter = parse_limiter(*limiter);
        cfg.limiter_name = to_lower(*limiter);
    }
    if (cfg.solver.recon.scheme == ReconstructionScheme::MUSCL && cfg.limiter_name.empty()) {
        throw std::runtime_error("reconstruction muscl requires a limiter");
    }

    const Grid grid(cfg.n_cells, cfg.x_min, cfg.x_max);
    const IdealGas gas(cfg.gamma);
    std::vector<Conserved> U0 = make_initial_condition(cfg.problem, grid, gas);

    std::filesystem::create_directories(cfg.output_dir);
    const std::filesystem::path dir(cfg.output_dir);
    write_solution_csv((dir / "solution_initial.csv").string(), U0, grid, gas);

    Solver solver(grid, gas, cfg.solver, std::move(U0));
    ResidualWriter residuals((dir / "residuals.csv").string());

    const auto start = std::chrono::steady_clock::now();
    const RunStats stats = solver.run([&](const StepInfo& info) { residuals.write(info); });
    const double wall =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    write_solution_csv((dir / "solution_final.csv").string(), solver.state(), grid, gas);

    std::optional<ErrorNorms> errors;
    if (cfg.problem.type == ProblemType::DensityWave) {
        errors = density_wave_errors(solver.state(), grid, cfg.problem.density_wave,
                                     stats.t_final_reached);
    }

    const bool muscl = cfg.solver.recon.scheme == ReconstructionScheme::MUSCL;
    SummaryInfo summary{cfg.problem_name,
                        cfg.n_cells,
                        grid.dx(),
                        cfg.flux_name,
                        cfg.reconstruction_name,
                        muscl && !cfg.limiter_name.empty() ? cfg.limiter_name : "none",
                        cfg.solver.cfl,
                        cfg.solver.t_final,
                        stats.n_steps,
                        wall,
                        errors};
    write_summary_csv((dir / "summary.csv").string(), summary);

    if (!quiet) {
        std::printf("problem=%s n_cells=%d n_steps=%d t=%.6f wall_time=%.3fs",
                    cfg.problem_name.c_str(), cfg.n_cells, stats.n_steps, stats.t_final_reached,
                    wall);
        if (errors) {
            std::printf(" L1_rho=%.6e L2_rho=%.6e Linf_rho=%.6e", errors->L1, errors->L2,
                        errors->Linf);
        }
        std::printf("\n");
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run_cli(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "euler1d: error: " << e.what() << std::endl;
        return 1;
    }
}
