#include <cstdlib>
#include <stdexcept>

#include "gryphon.h"

using namespace gryphon;

namespace {

void dumpFlux(const core::Input& input, const core::CosmicRays& cr) {
  const double units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;
  utils::OutputFile out(input.simname() + "_" + std::to_string(input.seed()) + ".txt");
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
    if (argc != 2) throw std::runtime_error("Usage: ./run seed");
    utils::Timer timer("timer for main");
    const auto seed = utils::parseSeed(argv[1]);

    auto in = core ::Input();
    in.set_seed(seed);
    in.set_simEmin(0.1 * cgs::TeV);
    in.set_simEmax(1e3 * cgs::TeV);
    in.set_simEsize(4 * 16);
    in.set_maxtime(100 * cgs::Myr);
    in.set_halosize(2 * cgs::kpc);
    in.set_simname("BASE_H2");
    in.set_spiralModel(SpiralModel::Steiman2010);
    in.print();

    RandomNumberGenerator rng(in.seed());

    auto galaxyModel = galaxy::makeGalaxy(in);
    galaxyModel->generate(rng);

    const auto& events = galaxyModel->get_events();

    LOGD << "event size : " << events.size();

    auto kernel = kernel::makeGreenKernel(in);
    auto injectionSpectra = injection::makeInjectionSpectra(in, events, rng);
    core::CosmicRays cr(in, kernel, std::move(injectionSpectra), events);
    cr.run();
    dumpFlux(in, cr);
  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
