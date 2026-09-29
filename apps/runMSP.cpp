#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "gryphon.h"
#include "gryphon/injection/MSP.h"
#include "gryphon/kernel/ContinuousLosses.h"

using namespace gryphon;

namespace {
struct Source {
  utils::Vector3d position;
  injection::MSPSpin spin;
};

void runMSP(const core::Input& input, const core::RunOutput& output) {
  RandomNumberGenerator rng(input.seed());
  auto galaxy = galaxy::makeGalaxy(input);
  injection::MSPInjection injection(input);
  kernel::ContinuousLosses solver(input);
  kernel::DiffusionLossesKernel transport(input);

  // An extant population, not a SN rate integrated for max_time. No source age
  // or random burst time is assigned: all sources inject steadily at current Lsd.
  std::vector<Source> sources;
  sources.reserve(input.mspSources());
  auto population = output.open("population", "x_sun [kpc] | y_sun [kpc] | z [kpc] | P [s] | B [G] | Pdot [s/s] | Lsd [erg/s]");
  population << std::setprecision(12);
  size_t outside = 0, nearby = 0;
  double totalPower = 0.;
  for (ulong i = 0; i < input.mspSources(); ++i) {
    Source s{galaxy->sample_position(rng), injection.sample(rng)};
    population << s.position.x/cgs::kpc << " " << s.position.y/cgs::kpc << " "
               << s.position.z/cgs::kpc << " " << s.spin.period/cgs::second << " "
               << s.spin.field/cgs::gauss << " " << s.spin.periodDot << " "
               << s.spin.luminosity/(cgs::erg/cgs::second) << "\n";
    totalPower += s.spin.luminosity;
    outside += std::abs(s.position.z) >= input.H();
    nearby += s.position.getModule() < cgs::kpc;
    sources.push_back(s);
  }
  std::vector<size_t> distanceOrder(sources.size());
  std::iota(distanceOrder.begin(), distanceOrder.end(), 0);
  std::sort(distanceOrder.begin(), distanceOrder.end(), [&](size_t i, size_t j) {
    return sources[i].position.getModule() < sources[j].position.getModule();
  });
  auto flux = output.open("flux", "E [GeV] | J_positron [GeV^-1 m^-2 s^-1 sr^-1]");
  auto stats = output.open("contributions", "E [GeV] | Neff | brightest_fraction | d50 [kpc] | d90 [kpc]");
  auto horizon = output.open("horizon", "E [GeV] | distance [kpc] | cumulative_fraction");
  auto checks = output.open("checks", "E [GeV] | max_relative_table_vs_direct (nearest 16) | relative_population_sum_change (double quadrature and distance grid)");
  flux << std::setprecision(12);
  stats << std::setprecision(12);
  horizon << std::setprecision(12);
  checks << std::setprecision(12);
  auto fineInput = input;
  fineInput.set_mspDistances(2*input.mspDistances()-1);
  fineInput.set_mspQuadrature(2*input.mspQuadrature());
  kernel::ContinuousLosses fine(fineInput);

  double largestChange = 0.;
  for (size_t e = 0; e < solver.energies().size(); ++e) {
    const double energy = solver.energies()[e];
    const double eGeV = energy/cgs::GeV;
    std::vector<double> values(sources.size());
    double sum = 0., squares = 0., largest = 0., refined = 0.;
    for (size_t i = 0; i < sources.size(); ++i) {
      const auto& s = sources[i];
      values[i] = s.spin.luminosity*solver.response(e, s.position);
      sum += values[i];
      squares += pow2(values[i]);
      largest = std::max(largest, values[i]);
      refined += s.spin.luminosity*fine.response(e, s.position);
    }
    const double change = sum > 0. ? std::abs(refined/sum-1.) : 0.;
    largestChange = std::max(largestChange, change);
    if (change > .005) throw std::runtime_error("MSP spectrum failed 0.5% convergence check");
    double directError = 0.;
    for (size_t j = 0; j < std::min<size_t>(16, sources.size()); ++j) {
      const auto& pos = sources[distanceOrder[j]].position;
      const double direct = solver.directResponse(energy, pos);
      if (direct > 0.) directError = std::max(directError, std::abs(solver.response(e, pos)/direct-1.));
    }
    if (directError > .005) throw std::runtime_error("MSP table/direct-kernel comparison exceeds 0.5%");
    checks << eGeV << " " << directError << " " << change << "\n";
    flux << eGeV << " " << sum*(cgs::GeV*cgs::m2*cgs::second*cgs::sr) << "\n";
    double cumulative = 0., d50 = 0., d90 = 0.;
    const bool diagnostic = std::abs(eGeV-10.) < 1e-8 || std::abs(eGeV-100.) < 1e-8 ||
                            std::abs(eGeV-300.) < 1e-8 || std::abs(eGeV-1000.) < 1e-8;
    for (size_t rank = 0; rank < distanceOrder.size(); ++rank) {
      const size_t i = distanceOrder[rank];
      cumulative += values[i];
      const double distance = sources[i].position.getModule()/cgs::kpc;
      if (sum > 0. && d50 == 0. && cumulative >= .5*sum) d50 = distance;
      if (sum > 0. && d90 == 0. && cumulative >= .9*sum) d90 = distance;
      if (diagnostic && (rank%50 == 0 || rank+1 == sources.size()))
        horizon << eGeV << " " << distance << " " << (sum > 0. ? cumulative/sum : 0.) << "\n";
    }
    stats << eGeV << " " << (squares > 0. ? sum*sum/squares : 0.) << " "
          << (sum > 0. ? largest/sum : 0.) << " " << d50 << " " << d90 << "\n";
  }
  auto summary = output.open("summary", "quantity = value (CGS unless indicated)");
  summary << std::setprecision(12)
          << "calculation = stationary_MSP_positrons\n"
          << "sources = " << sources.size() << "\n"
          << "outside_halo = " << outside << "\n"
          << "sources_within_1kpc = " << nearby << "\n"
          << "total_spin_down_erg_s = " << totalPower/(cgs::erg/cgs::second) << "\n"
          << "b_at_10GeV_GeV_Myr = " << transport.b(10.*cgs::GeV)/(cgs::GeV/cgs::Myr) << "\n"
          << "image_order = " << solver.imageOrder() << "\n"
          << "max_relative_refinement = " << largestChange << "\n";
  LOGI << "MSPs: " << sources.size() << "; total spin-down = " << totalPower
       << " erg/s; maximum refinement change = " << largestChange;
}
}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "runMSP");
    if (!args) return EXIT_SUCCESS;
    const auto input = utils::makeInput(*args);
    if (input.injectionModel() != InjectionModel::MSP)
      throw std::invalid_argument("runMSP requires injectionmodel=MSP");
    if (input.mspQuadrature() > 2048 || input.mspDistances() > 50000)
      throw std::invalid_argument("runMSP needs quadrature <= 2048 and distances <= 50000 to run the doubled-resolution check");
    // Do not silently overwrite an existing realization.
    const auto dir = std::filesystem::path(args->outdir)/input.simname();
    std::ostringstream seed;
    seed << std::setfill('0') << std::setw(6) << input.seed();
    if (std::filesystem::exists(dir/("flux_"+seed.str()+".txt")))
      throw std::runtime_error("MSP realization already exists; choose another seed or --outdir");
    input.print();
    utils::Timer timer("runMSP");
    const core::RunOutput output(input, args->outdir);
    // Keep the parameters for each seed even if another run replaces params.ini.
    input.write_params_file(output.path("params").string());
    runMSP(input, output);
  } catch (const std::exception& error) {
    LOGE << "!Fatal Error: " << error.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
