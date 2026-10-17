#include "euler1d/io.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <map>
#include <set>
#include <stdexcept>

namespace euler1d {

namespace {

struct Entry {
    std::string value;
    int line;
};
using Table = std::map<std::string, Entry>;

const std::set<std::string> kGeneralKeys = {
    "problem", "n_cells", "x_min", "x_max", "t_final", "cfl", "flux", "reconstruction",
    "bc_left", "bc_right", "output_dir", "gamma", "limiter", "max_steps"};
const std::set<std::string> kRiemannKeys = {"x0",    "rho_L", "u_L", "p_L",
                                            "rho_R", "u_R",   "p_R"};
const std::set<std::string> kDensityWaveKeys = {"rho0", "amplitude", "u0", "p0", "wavenumber"};

[[noreturn]] void fail(int line, const std::string& msg) {
    throw std::runtime_error("config line " + std::to_string(line) + ": " + msg);
}

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

const Entry& require(const Table& t, const std::string& key) {
    const auto it = t.find(key);
    if (it == t.end()) throw std::runtime_error("config: missing required key '" + key + "'");
    return it->second;
}

double to_double(const Entry& e, const std::string& key) {
    try {
        std::size_t pos = 0;
        const double v = std::stod(e.value, &pos);
        if (pos != e.value.size()) throw std::invalid_argument("trailing characters");
        return v;
    } catch (const std::exception&) {
        fail(e.line, "cannot parse '" + e.value + "' as a number for key '" + key + "'");
    }
}

int to_int(const Entry& e, const std::string& key) {
    try {
        std::size_t pos = 0;
        const int v = std::stoi(e.value, &pos);
        if (pos != e.value.size()) throw std::invalid_argument("trailing characters");
        return v;
    } catch (const std::exception&) {
        fail(e.line, "cannot parse '" + e.value + "' as an integer for key '" + key + "'");
    }
}

double get_double(const Table& t, const std::string& key) {
    return to_double(require(t, key), key);
}

int get_int(const Table& t, const std::string& key) { return to_int(require(t, key), key); }

template <typename Parser>
auto get_enum(const Table& t, const std::string& key, Parser parse) {
    const Entry& e = require(t, key);
    try {
        return parse(e.value);
    } catch (const std::invalid_argument& ex) {
        fail(e.line, ex.what());
    }
}

void reject_keys(const Table& t, const std::set<std::string>& keys, const std::string& problem) {
    for (const auto& key : keys) {
        const auto it = t.find(key);
        if (it != t.end()) {
            fail(it->second.line, "key '" + key + "' is not valid for problem '" + problem + "'");
        }
    }
}

std::string csv_double(double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.16e", v);
    return buf;
}

void ensure_parent_dir(const std::string& path) {
    const auto parent = std::filesystem::path(path).parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
}

std::ofstream open_output(const std::string& path) {
    ensure_parent_dir(path);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open '" + path + "' for writing");
    return out;
}

}  // namespace

Config parse_config(std::istream& in) {
    Table table;
    std::string raw;
    int line_no = 0;
    while (std::getline(in, raw)) {
        ++line_no;
        const auto hash = raw.find('#');
        if (hash != std::string::npos) raw.erase(hash);
        const std::string line = trim(raw);
        if (line.empty()) continue;

        const auto eq = line.find('=');
        if (eq == std::string::npos) fail(line_no, "expected 'key = value'");
        const std::string key = trim(line.substr(0, eq));
        const std::string value = trim(line.substr(eq + 1));
        if (key.empty() || value.empty()) fail(line_no, "expected 'key = value'");

        if (!kGeneralKeys.count(key) && !kRiemannKeys.count(key) &&
            !kDensityWaveKeys.count(key)) {
            fail(line_no, "unknown key '" + key + "'");
        }
        if (table.count(key)) fail(line_no, "duplicate key '" + key + "'");
        table[key] = {value, line_no};
    }

    Config cfg{};
    const std::string problem_value = require(table, "problem").value;
    cfg.problem.type = get_enum(table, "problem", parse_problem_type);
    cfg.problem_name = to_lower(problem_value);

    if (cfg.problem.type == ProblemType::Riemann) {
        reject_keys(table, kDensityWaveKeys, cfg.problem_name);
        RiemannParams& r = cfg.problem.riemann;
        r = {get_double(table, "x0"),    get_double(table, "rho_L"), get_double(table, "u_L"),
             get_double(table, "p_L"),   get_double(table, "rho_R"), get_double(table, "u_R"),
             get_double(table, "p_R")};
    } else {
        reject_keys(table, kRiemannKeys, cfg.problem_name);
        cfg.problem.density_wave = {get_double(table, "rho0"), get_double(table, "amplitude"),
                                    get_double(table, "u0"), get_double(table, "p0"),
                                    get_int(table, "wavenumber")};
    }

    cfg.n_cells = get_int(table, "n_cells");
    cfg.x_min = get_double(table, "x_min");
    cfg.x_max = get_double(table, "x_max");
    cfg.solver.t_final = get_double(table, "t_final");
    cfg.solver.cfl = get_double(table, "cfl");
    cfg.output_dir = require(table, "output_dir").value;

    cfg.solver.flux = get_enum(table, "flux", parse_flux_scheme);
    cfg.flux_name = to_lower(require(table, "flux").value);
    cfg.solver.recon.scheme = get_enum(table, "reconstruction", parse_reconstruction);
    cfg.reconstruction_name = to_lower(require(table, "reconstruction").value);
    cfg.solver.left = get_enum(table, "bc_left", parse_boundary_type);
    cfg.solver.right = get_enum(table, "bc_right", parse_boundary_type);

    cfg.solver.recon.limiter = Limiter::Minmod;
    if (table.count("limiter")) {
        cfg.solver.recon.limiter = get_enum(table, "limiter", parse_limiter);
        cfg.limiter_name = to_lower(table.at("limiter").value);
    } else if (cfg.solver.recon.scheme == ReconstructionScheme::MUSCL) {
        throw std::runtime_error("config: missing required key 'limiter' (reconstruction = muscl)");
    }

    cfg.gamma = table.count("gamma") ? get_double(table, "gamma") : 1.4;
    cfg.solver.max_steps = table.count("max_steps") ? get_int(table, "max_steps") : 1000000;
    return cfg;
}

Config load_config(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open config file '" + path + "'");
    return parse_config(in);
}

void write_solution_csv(const std::string& path, const std::vector<Conserved>& U,
                        const Grid& grid, const IdealGas& gas) {
    std::ofstream out = open_output(path);
    out << "x,rho,u,p,e\n";
    for (int i = grid.begin(); i < grid.end(); ++i) {
        const Primitive W = gas.to_primitive(U[i]);
        const double e = W.p / ((gas.gamma() - 1.0) * W.rho);
        out << csv_double(grid.x_center(i)) << ',' << csv_double(W.rho) << ','
            << csv_double(W.u) << ',' << csv_double(W.p) << ',' << csv_double(e) << '\n';
    }
}

ResidualWriter::ResidualWriter(const std::string& path) : out_(open_output(path)) {
    out_ << "step,t,dt,res_rho,res_mom,res_E\n";
}

void ResidualWriter::write(const StepInfo& info) {
    out_ << info.step << ',' << csv_double(info.t) << ',' << csv_double(info.dt) << ','
         << csv_double(info.residual_L2[0]) << ',' << csv_double(info.residual_L2[1]) << ','
         << csv_double(info.residual_L2[2]) << '\n';
}

void write_summary_csv(const std::string& path, const SummaryInfo& info) {
    std::ofstream out = open_output(path);
    out << "problem,n_cells,dx,flux,reconstruction,limiter,cfl,t_final,n_steps,wall_time_s,"
           "L1_rho,L2_rho,Linf_rho\n";
    out << info.problem << ',' << info.n_cells << ',' << csv_double(info.dx) << ','
        << info.flux << ',' << info.reconstruction << ',' << info.limiter << ','
        << csv_double(info.cfl) << ',' << csv_double(info.t_final) << ',' << info.n_steps << ','
        << csv_double(info.wall_time_s) << ',';
    if (info.rho_errors) {
        out << csv_double(info.rho_errors->L1) << ',' << csv_double(info.rho_errors->L2) << ','
            << csv_double(info.rho_errors->Linf);
    } else {
        out << ",,";
    }
    out << '\n';
}

}  // namespace euler1d
