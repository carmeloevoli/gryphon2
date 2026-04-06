#include <array>
#include <cstdlib>
#include <stdexcept>

#include "gryphon.h"

using namespace gryphon;

namespace {

struct ModelScenario {
  const char* sim_name;
  SpiralModel spiral_model;
  double halo_size;
  double sn_rate;
  bool vary_energy;
  bool vary_slope;
};

core::Input makeInput(unsigned long int seed, const ModelScenario& scenario) {
  auto input = core::Input();
  input.set_seed(seed);
  input.set_simname(scenario.sim_name);
  input.set_spiralModel(scenario.spiral_model);
  input.set_injectionModel(InjectionModel::GalacticRandom);
  input.set_injEmax(1e12 * cgs::GeV);
  input.set_simEmin(1e2 * cgs::GeV);
  input.set_simEmax(1e6 * cgs::GeV);
  input.set_simEsize(16 * 4);
  input.set_maxtime(100. * cgs::Myr);

  const auto base_halo_size = 4. * cgs::kpc;
  input.set_halosize(scenario.halo_size * base_halo_size);

  const auto base_sn_rate = 1. / 50. / cgs::year;
  input.set_rate(scenario.sn_rate * base_sn_rate);

  if (scenario.vary_energy) {
    input.enable_varyenergy();
  }
  if (scenario.vary_slope) {
    input.enable_varyslope();
  }

  return input;
}

constexpr std::array<ModelScenario, 6> kModelScenarios{{
    {"test_fixed2", SpiralModel::Steiman2010, 0.5, 1.0, false, false},
    {"test_fixed4", SpiralModel::Steiman2010, 1.0, 1.0, false, false},
    {"test_fixed8", SpiralModel::Steiman2010, 2.0, 1.0, false, false},
    {"test_raresn", SpiralModel::Steiman2010, 1.0, 0.1, false, false},
    {"test_varyesn", SpiralModel::Steiman2010, 1.0, 1.0, true, false},
    {"test_varyinj", SpiralModel::Steiman2010, 1.0, 1.0, false, true},
}};

void runAndDump(const core::Input& in) {
  RandomNumberGenerator rng(in.seed());
  auto galaxyModel = galaxy::makeGalaxy(in);
  galaxyModel->generate(rng, false);
  auto events = galaxyModel->get_events();

  auto kernel = kernel::makeGreenKernel(in);
  auto injectionSpectra = injection::makeInjectionSpectra(in, events, rng);
  core::CosmicRays cr(in, kernel, std::move(injectionSpectra), events);
  cr.run();

  const double units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;
  utils::OutputFile out(in.simname() + "_" + std::to_string(in.seed()) + ".txt");
  out << "# E [GeV] - I [GeV-1 m-2 sec-1 sr-1]\n";
  out << std::scientific;
  const auto& E = cr.get_energyAxis();
  const auto& I = cr.get_flux();
  for (size_t i = 0; i < E.size(); ++i) {
    out << E[i] / cgs::GeV << "\t";
    out << I[i] / units << "\t\n";
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    if (argc != 2) throw std::runtime_error("Usage: ./runModels seed");
    utils::Timer timer("timer for main");
    const auto seed = utils::parseSeed(argv[1]);

    for (const auto& scenario : kModelScenarios) {
      runAndDump(makeInput(seed, scenario));
    }

  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
