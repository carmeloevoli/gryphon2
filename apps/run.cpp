#include "gryphon.h"

using namespace gryphon;

namespace {

void dumpFlux(const core::Input& input, const core::CosmicRays& cr) {
  const double units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;
  utils::OutputFile out(input.simname + "_" + std::to_string(input.seed) + ".txt");
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
    if (argc != 2) throw std::runtime_error("Usage: ./run params.ini");
    utils::Timer timer("timer for main");

    // TODO check if output dir exists

    auto in = core ::Input();
    in.set_seed(atoi(argv[1]));
    in.set_simEmin(0.1 * cgs::TeV);
    in.set_simEmax(1e3 * cgs::TeV);
    in.set_simEsize(4 * 16);
    in.set_maxtime(100 * cgs::Myr);
    in.set_halosize(2 * cgs::kpc);
    in.set_simname("base_H2");
    in.set_spiralModel(SpiralModel::Steiman2010);
    in.print();

    RandomNumberGenerator rng = utils::RNG<double>(in.seed);

    std::shared_ptr<galaxy::Galaxy> galaxy;
    switch (in.spiralModel) {
      case SpiralModel::Uniform:
        galaxy = std::make_shared<galaxy::GalaxyUniform>(in);
        break;
      case SpiralModel::Jelly:
        galaxy = std::make_shared<galaxy::GalaxyJelly>(in);
        break;
      case SpiralModel::Steiman2010:
        galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
        break;
      default:
        throw std::invalid_argument("Spiral model not implemented yet");
    }
    galaxy->generate(rng);

    auto events = galaxy->get_events();

    LOGD << "event size : " << events.size();

    auto kernel = kernel::makeGreenKernel(in);
    auto injectionSpectrum = injection::makeInjectionSpectrum(in, rng);
    core::CosmicRays cr(in, kernel, injectionSpectrum, events);
    cr.run();
    dumpFlux(in, cr);
  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
  }
  return EXIT_SUCCESS;
}