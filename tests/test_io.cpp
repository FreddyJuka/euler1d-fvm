#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "euler1d/io.hpp"

using namespace euler1d;

namespace {

const std::string kGeneral =
    "n_cells = 100\nx_min = 0.0\nx_max = 1.0\nt_final = 0.2\ncfl = 0.8\n"
    "flux = HLLC\nreconstruction = muscl\nlimiter = minmod\n"
    "bc_left = transmissive\nbc_right = transmissive\noutput_dir = out\n";
const std::string kRiemann =
    "problem = riemann\nx0 = 0.5\nrho_L = 1.0\nu_L = 0.0\np_L = 1.0\n"
    "rho_R = 0.125\nu_R = 0.0\np_R = 0.1\n";
const std::string kWave =
    "problem = density_wave\nrho0 = 1.0\namplitude = 0.2\nu0 = 1.0\np0 = 1.0\nwavenumber = 1\n";

Config parse(const std::string& text) {
    std::istringstream in(text);
    return parse_config(in);
}

std::filesystem::path temp_dir(const std::string& name) {
    const auto dir = std::filesystem::temp_directory_path() / ("euler1d_test_" + name);
    std::filesystem::remove_all(dir);
    return dir;
}

int count_lines(const std::filesystem::path& p, std::string* header = nullptr) {
    std::ifstream in(p);
    std::string line;
    int n = 0;
    while (std::getline(in, line)) {
        if (n == 0 && header) *header = line;
        ++n;
    }
    return n;
}

}  // namespace

TEST(Config, ValidRiemann) {
    const Config c = parse("# comment\n\n" + kRiemann + kGeneral + "gamma = 1.67  # inline\nmax_steps = 50\n");
    EXPECT_EQ(c.problem.type, ProblemType::Riemann);
    EXPECT_EQ(c.problem.riemann.x0, 0.5);
    EXPECT_EQ(c.problem.riemann.rho_R, 0.125);
    EXPECT_EQ(c.problem.riemann.p_R, 0.1);
    EXPECT_EQ(c.n_cells, 100);
    EXPECT_EQ(c.x_max, 1.0);
    EXPECT_EQ(c.gamma, 1.67);
    EXPECT_EQ(c.solver.max_steps, 50);
    EXPECT_EQ(c.solver.cfl, 0.8);
    EXPECT_EQ(c.solver.t_final, 0.2);
    EXPECT_EQ(c.solver.flux, FluxScheme::HLLC);
    EXPECT_EQ(c.solver.recon.scheme, ReconstructionScheme::MUSCL);
    EXPECT_EQ(c.solver.recon.limiter, Limiter::Minmod);
    EXPECT_EQ(c.solver.left, BoundaryType::Transmissive);
    EXPECT_EQ(c.output_dir, "out");
}

TEST(Config, ValidDensityWaveDefaults) {
    const Config c = parse(kWave + kGeneral);
    EXPECT_EQ(c.problem.type, ProblemType::DensityWave);
    EXPECT_EQ(c.problem.density_wave.amplitude, 0.2);
    EXPECT_EQ(c.problem.density_wave.wavenumber, 1);
    EXPECT_EQ(c.gamma, 1.4);
    EXPECT_EQ(c.solver.max_steps, 1000000);
}

TEST(Config, Errors) {
    EXPECT_THROW(parse(kRiemann + kGeneral + "bogus = 1\n"), std::runtime_error);
    EXPECT_THROW(parse(kRiemann + kGeneral + "cfl = 0.5\n"), std::runtime_error);
    EXPECT_THROW(parse(kRiemann + "n_cells = 10\n"), std::runtime_error);
    EXPECT_THROW(parse(kRiemann + kGeneral + "gamma = abc\n"), std::runtime_error);
    EXPECT_THROW(parse(kRiemann + kGeneral + "this is not a pair\n"), std::runtime_error);
    EXPECT_THROW(parse(kWave + kRiemann.substr(kRiemann.find('\n') + 1) + kGeneral),
                 std::runtime_error);

    std::string no_limiter = kGeneral;
    no_limiter.erase(no_limiter.find("limiter = minmod\n"), std::string("limiter = minmod\n").size());
    EXPECT_THROW(parse(kRiemann + no_limiter), std::runtime_error);
}

TEST(Config, RejectsMalformedValues) {
    const auto with = [](const std::string& from, const std::string& to) {
        std::string text = kRiemann + kGeneral;
        const auto pos = text.find(from);
        EXPECT_NE(pos, std::string::npos);
        text.replace(pos, from.size(), to);
        return text;
    };
    EXPECT_THROW(parse(with("n_cells = 100", "n_cells = 1e3")), std::runtime_error);
    EXPECT_THROW(parse(with("problem = riemann", "problem = bogus_problem")), std::runtime_error);
    EXPECT_THROW(parse(with("flux = HLLC", "flux = bogus_flux")), std::runtime_error);
    EXPECT_THROW(parse(kRiemann + kGeneral + "x0 = 0.4\n"), std::runtime_error);  // duplicate key
}

TEST(Config, ErrorReportsLineNumber) {
    try {
        parse("problem = riemann\n\nbogus = 1\n");
        FAIL() << "expected exception";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("line 3"), std::string::npos);
    }
}

TEST(Config, LoadMissingFileThrows) {
    EXPECT_THROW(load_config("/nonexistent/path.cfg"), std::runtime_error);
}

TEST(Csv, Writers) {
    const auto dir = temp_dir("csv");
    const Grid grid(8, 0.0, 1.0);
    const IdealGas gas(1.4);
    std::vector<Conserved> U(grid.size(), gas.to_conserved({1.0, 0.5, 1.0}));

    std::string header;
    write_solution_csv((dir / "sol.csv").string(), U, grid, gas);
    EXPECT_EQ(count_lines(dir / "sol.csv", &header), 1 + 8);
    EXPECT_EQ(header, "x,rho,u,p,e");

    {
        ResidualWriter w((dir / "res.csv").string());
        w.write({1, 0.1, 0.1, {1.0, 2.0, 3.0}});
        w.write({2, 0.2, 0.1, {1.0, 2.0, 3.0}});
    }
    EXPECT_EQ(count_lines(dir / "res.csv", &header), 1 + 2);
    EXPECT_EQ(header, "step,t,dt,res_rho,res_mom,res_E");

    SummaryInfo s{"density_wave", 8, 0.125, "hllc", "muscl", "vanleer", 0.8, 1.0, 1.0, 10, 0.01,
                  ErrorNorms{1.0, 2.0, 3.0}};
    write_summary_csv((dir / "summary.csv").string(), s);
    EXPECT_EQ(count_lines(dir / "summary.csv", &header), 2);
    EXPECT_EQ(header,
              "problem,n_cells,dx,flux,reconstruction,limiter,cfl,t_final,t_reached,n_steps,wall_time_s,"
              "L1_rho,L2_rho,Linf_rho");
    std::filesystem::remove_all(dir);
}

TEST(Config, RejectsNonFiniteNumbers) {
    for (const char* bad : {"nan", "inf"}) {
        const std::string text = "problem = riemann\nx0 = " + std::string(bad) +
                                 "\nrho_L = 1.0\nu_L = 0.0\np_L = 1.0\n"
                                 "rho_R = 0.125\nu_R = 0.0\np_R = 0.1\n" + kGeneral;
        EXPECT_THROW(parse(text), std::runtime_error) << bad;
    }
    std::string text = kRiemann + kGeneral;
    text.replace(text.find("t_final = 0.2"), 13, "t_final = inf");
    EXPECT_THROW(parse(text), std::runtime_error);
}
