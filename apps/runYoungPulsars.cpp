#include <cstdlib>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>

#include "gryphon.h"

using namespace gryphon;

namespace {

core::Input makeInput(unsigned long int seed) {
  auto input = core::Input();
  input.set_seed(seed);
  input.set_simname("youngpulsars_H4_100msec");
  input.set_pid(core::H);

  input.set_spiralModel(SpiralModel::Steiman2010);
  input.set_transportModel(TransportModel::PureDiffusion);
  input.set_injectionModel(InjectionModel::YoungPulsars);

  input.set_simEmin(1e4 * cgs::GeV);
  input.set_simEmax(1e8 * cgs::GeV);
  input.set_simEsize(16 * 4);
  input.set_maxtime(100. * cgs::Myr);

  input.set_galaxyRadius(20. * cgs::kpc);
  input.set_sunRadius(8.3 * cgs::kpc);
  input.set_halosize(4. * cgs::kpc);
  input.set_rate(1. / 60. / cgs::year);

  input.set_efficiency(1.);
  input.set_youngPulsarsP0(100. * cgs::msec);
  input.set_youngPulsarsSigmaP0(10. * cgs::msec);
  input.set_youngPulsarsRandomInitialPeriod(true);
  input.set_youngPulsarsB0(2.5e12 * cgs::gauss);
  input.set_youngPulsarsSigmaLog10B(0.5);
  input.set_youngPulsarsRandomMagneticField(true);

  return input;
}

void dumpFlux(const core::Input& input, const core::CosmicRays& cr) {
  const double flux_units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;
  utils::OutputFile out(input.simname() + "_" + std::to_string(input.seed()) + ".txt");
  out << "# E [GeV] - I [GeV^-1 m^-2 s^-1 sr^-1]\n";
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

void runYoungPulsars(const core::Input& input) {
  RandomNumberGenerator rng(input.seed());

  auto galaxyModel = galaxy::makeGalaxy(input);
  galaxyModel->generate(rng, false);
  const auto& events = galaxyModel->get_events();
  LOGD << "event size : " << events.size();

  auto kernel = kernel::makeGreenKernel(input);
  auto injectionSpectra = injection::makeInjectionSpectra(input, events, rng);
  core::CosmicRays cr(input, kernel, std::move(injectionSpectra), events);
  cr.run();

  dumpFlux(input, cr);
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    if (argc != 2) throw std::runtime_error("Usage: ./runYoungPulsars seed");
    utils::Timer timer("timer for main");

    const auto input = makeInput(utils::parseSeed(argv[1]));
    input.print();
    runYoungPulsars(input);

  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
