#include <cstdlib>
#include <iomanip>
#include <stdexcept>
#include <string>

#include "gryphon.h"

using namespace gryphon;

namespace {

core::Input makeInput(unsigned long int seed) {
  auto input = core::Input();
  input.set_seed(seed);
  input.set_simname("randomemax_H4_5PeV_sigma0p4");
  input.set_pid(core::H);

  input.set_spiralModel(SpiralModel::Steiman2010);
  input.set_transportModel(TransportModel::PureDiffusion);
  input.set_injectionModel(InjectionModel::RandomEmax);

  input.set_simEmin(1e2 * cgs::GeV);
  input.set_simEmax(1e8 * cgs::GeV);
  input.set_simEsize(16 * 4);
  input.set_maxtime(100. * cgs::Myr);

  input.set_galaxyRadius(20. * cgs::kpc);
  input.set_sunRadius(8.3 * cgs::kpc);
  input.set_halosize(4. * cgs::kpc);
  input.set_rate(1. / 50. / cgs::year);

  input.set_efficiency(0.1);
  input.set_injSlope(2.32);
  input.set_injEmax(5. * cgs::PeV);
  input.set_injEmaxSigmaDex(0.4);
  input.set_injEmaxMin(100. * cgs::TeV);
  input.set_injEmaxMax(100. * cgs::PeV);

  return input;
}

std::string outputStem(const core::Input& input) {
  return input.simname() + "_" + std::to_string(input.seed());
}

void dumpFlux(const core::Input& input, const core::CosmicRays& cr) {
  const double flux_units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;
  utils::OutputFile out(outputStem(input) + ".txt");
  out << "# E [GeV] - I [GeV^-1 m^-2 s^-1 sr^-1]\n";
  out << std::scientific << std::setprecision(6);

  const auto& E = cr.get_energyAxis();
  const auto& I = cr.get_flux();
  for (size_t i = 0; i < E.size(); ++i) {
    out << E[i] / cgs::GeV << "\t";
    out << I[i] / flux_units << "\n";
  }
}

void dumpSourceCatalog(const core::Input& input, const core::Events& events,
                       const injection::InjectionSpectra& spectra) {
  if (events.size() != spectra.size()) {
    throw std::runtime_error("RandomEmax source catalog requires one spectrum per event");
  }

  const double invMyr = 1. / cgs::Myr;
  const double invKpc = 1. / cgs::kpc;
  const double invTeV = 1. / cgs::TeV;
  const double invErg = 1. / cgs::erg;

  utils::OutputFile out(outputStem(input) + "_sources.txt");
  out << "# age [Myr]\tx [kpc]\ty [kpc]\tz [kpc]\tr [kpc]\tEmax [TeV]\tcrEnergy [erg]\n";
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

void runRandomEmaxPopulation(const core::Input& input) {
  RandomNumberGenerator rng(input.seed());

  auto galaxyModel = galaxy::makeGalaxy(input);
  galaxyModel->generate(rng, false);
  const auto& events = galaxyModel->get_events();
  LOGD << "event size : " << events.size();

  auto kernel = kernel::makeGreenKernel(input);
  auto injectionSpectra = injection::makeInjectionSpectra(input, events, rng);

  dumpSourceCatalog(input, events, injectionSpectra);

  core::CosmicRays cr(input, kernel, std::move(injectionSpectra), events);
  cr.run();

  dumpFlux(input, cr);
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    if (argc != 2) throw std::runtime_error("Usage: ./runRandomEmax seed");
    utils::Timer timer("timer for main");

    const auto input = makeInput(utils::parseSeed(argv[1]));
    input.print();
    runRandomEmaxPopulation(input);

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
