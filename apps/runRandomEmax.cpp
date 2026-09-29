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
    out << E[i] / cgs::GeV << "\t";
    out << I[i] / flux_units << "\n";
  }
}

void dumpSourceCatalog(const core::RunOutput& run, const core::Events& events,
                       const injection::InjectionSpectra& spectra) {
  if (events.size() != spectra.size()) {
    throw std::runtime_error("RandomEmax source catalog requires one spectrum per event");
  }

  const double invMyr = 1. / cgs::Myr;
  const double invKpc = 1. / cgs::kpc;
  const double invTeV = 1. / cgs::TeV;
  const double invErg = 1. / cgs::erg;

  auto out = run.open("sources",
                      "age [Myr] | x [kpc] | y [kpc] | z [kpc] | r [kpc] | "
                      "Emax [TeV] | crEnergy [erg]");
  out << std::scientific << std::setprecision(6);

  for (size_t i = 0; i < events.size(); ++i) {
    const auto& event = events[i];
    if (!event) continue;

    const auto* spectrum = dynamic_cast<const injection::RandomEmaxSpectrum*>(spectra[i].get());
    if (spectrum == nullptr) {
      throw std::runtime_error("Encountered a non-RandomEmax spectrum in runRandomEmax");
    }

    out << event->age * invMyr << "\t";
    out << event->pos.getX() * invKpc << "\t";
    out << event->pos.getY() * invKpc << "\t";
    out << event->pos.getZ() * invKpc << "\t";
    out << event->pos.getModule() * invKpc << "\t";
    out << spectrum->Emax * invTeV << "\t";
    out << spectrum->crEnergy * invErg << "\n";
  }
}

void runRandomEmaxPopulation(const core::Input& input, const core::RunOutput& run) {
  RandomNumberGenerator rng(input.seed());

  auto galaxyModel = galaxy::makeGalaxy(input);
  galaxyModel->generate(rng, false);
  const auto& events = galaxyModel->get_events();
  LOGD << "event size : " << events.size();

  auto kernel = kernel::makeGreenKernel(input);
  auto injectionSpectra = injection::makeInjectionSpectra(input, events, rng);

  dumpSourceCatalog(run, events, injectionSpectra);

  core::CosmicRays cr(input, kernel, std::move(injectionSpectra), events);
  cr.run();

  dumpFlux(run, cr);
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "runRandomEmax");
    if (!args) return EXIT_SUCCESS;

    utils::Timer timer("timer for main");
    const auto input = utils::makeInput(*args);
    input.print();

    const core::RunOutput run(input, args->outdir);
    runRandomEmaxPopulation(input, run);

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
