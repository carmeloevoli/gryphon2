#include <array>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <string_view>

#include "gryphon.h"

using namespace gryphon;

namespace {

struct ModelScenario {
  const char* sim_name;
  SpiralModel spiral_model;
  double halo_size;
  double sn_rate;
  double cr_energy_scale;
  bool tail_only;
  bool vary_energy;
  bool vary_slope;
};

core::Input makeInput(unsigned long int seed, const ModelScenario& scenario) {
  auto input = core::Input();
  input.set_seed(seed);
  input.set_simname(scenario.sim_name);
  input.set_spiralModel(scenario.spiral_model);
  input.set_injectionModel(InjectionModel::GalacticRandom);
  input.set_injEmax(0.);  // Emax -> infinity; indices conditioned on gamma > 2.
  input.set_injSlopeSigma(0.30);  // Width of the parent Gaussian before truncation.
  input.set_simEmin((scenario.tail_only ? 1e3 : 1e2) * cgs::GeV);
  input.set_simEmax(1e6 * cgs::GeV);
  input.set_simEsize(scenario.tail_only ? 48 : 16 * 4);
  input.set_maxtime((scenario.tail_only ? 30. : 100.) * cgs::Myr);

  const auto base_halo_size = 4. * cgs::kpc;
  input.set_halosize(scenario.halo_size * base_halo_size);

  const auto base_sn_rate = 1. / 50. / cgs::year;
  input.set_rate(scenario.sn_rate * base_sn_rate);
  input.set_efficiency(0.1 * scenario.cr_energy_scale);

  if (scenario.vary_energy) {
    input.enable_varyenergy();
  }
  if (scenario.vary_slope) {
    input.enable_varyslope();
  }

  return input;
}

constexpr std::array<ModelScenario, 6> kModelScenarios{{
    {"test_fixed2", SpiralModel::Xie2024, 0.5, 1.0, 1.0, false, false, false},
    {"test_fixed4", SpiralModel::Xie2024, 1.0, 1.0, 1.0, false, false, false},
    {"test_fixed8", SpiralModel::Xie2024, 2.0, 1.0, 1.0, false, false, false},
    // Ten per cent of remnants supply the high-energy population.  Their CR
    // energy is increased by ten so that the time-averaged CR power is fixed.
    // This tail-only run starts at 1 TeV and spans 30 Myr, more than six
    // escape times at its lowest energy for H=4 kpc.
    {"test_raresn", SpiralModel::Xie2024, 1.0, 0.1, 10.0, true, false, false},
    {"test_varyesn", SpiralModel::Xie2024, 1.0, 1.0, 1.0, false, true, false},
    {"test_varyinj", SpiralModel::Xie2024, 1.0, 1.0, 1.0, false, false, true},
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
    if (argc < 2 || argc > 3) {
      throw std::runtime_error("Usage: ./runModels seed [scenario]");
    }
    utils::Timer timer("timer for main");
    const auto seed = utils::parseSeed(argv[1]);
    const auto requested_scenario =
        (argc == 3) ? std::string_view(argv[2]) : std::string_view{};
    bool scenario_found = requested_scenario.empty();

    for (const auto& scenario : kModelScenarios) {
      if (!requested_scenario.empty() && requested_scenario != scenario.sim_name) continue;
      scenario_found = true;
      runAndDump(makeInput(seed, scenario));
    }
    if (!scenario_found) {
      throw std::runtime_error("Unknown model scenario: " + std::string(requested_scenario));
    }

  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
