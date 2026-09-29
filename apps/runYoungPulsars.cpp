#include <cstdlib>
#include <iomanip>
#include <stdexcept>

#include "gryphon.h"

using namespace gryphon;

namespace {

void dumpFlux(const core::RunOutput& run, const core::CosmicRays& cr) {
  const double flux_units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;
  auto out = run.open("flux", "E [GeV] | I [GeV^-1 m^-2 s^-1 sr^-1]");
  out << std::scientific << std::setprecision(6);

  const auto& E = cr.get_energyAxis();
  const auto& I = cr.get_flux();
  for (size_t i = 0; i < E.size(); ++i) {
    const double energy_gev = E[i] / cgs::GeV;
    const double flux = I[i] / flux_units;
    out << energy_gev << "\t";
    out << flux << "\t";
    out << "\n";
  }
}

void runYoungPulsars(const core::Input& input, const core::RunOutput& run) {
  RandomNumberGenerator rng(input.seed());

  auto galaxyModel = galaxy::makeGalaxy(input);
  galaxyModel->generate(rng, false);
  const auto& events = galaxyModel->get_events();
  LOGD << "event size : " << events.size();

  auto kernel = kernel::makeGreenKernel(input);
  auto injectionSpectra = injection::makeInjectionSpectra(input, events, rng);
  core::CosmicRays cr(input, kernel, std::move(injectionSpectra), events);
  cr.run();

  dumpFlux(run, cr);
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "runYoungPulsars");
    if (!args) return EXIT_SUCCESS;

    utils::Timer timer("timer for main");
    const auto input = utils::makeInput(*args);
    input.print();

    const core::RunOutput run(input, args->outdir);
    runYoungPulsars(input, run);

  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
