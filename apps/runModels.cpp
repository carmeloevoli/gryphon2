#include "gryphon.h"

using namespace gryphon;

namespace {

void runAndDump(const core::Input& in, const core::Events& events, RandomNumberGenerator& rng) {
  auto kernel = kernel::makeGreenKernel(in);
  auto injectionSpectrum = injection::makeInjectionSpectrum(in, rng);
  core::CosmicRays cr(in, kernel, injectionSpectrum, events);
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

void runConfigured(core::Input& in) {
  RandomNumberGenerator rng(in.seed());
  auto galaxyModel = galaxy::makeGalaxy(in);
  galaxyModel->generate(rng, false);
  runAndDump(in, galaxyModel->get_events(), rng);
}

}  // namespace

void runUniform(unsigned long int seed, std::string simname) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.set_simname(simname);
  in.set_spiralModel(SpiralModel::Uniform);
  // in.print();

  runConfigured(in);
}

void runFixed(unsigned long int seed, double haloSize, std::string simName) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.set_halosize(haloSize);
  in.set_simname(simName);
  in.set_spiralModel(SpiralModel::Steiman2010);
  // in.print();

  runConfigured(in);
}

void runVaryEnergy(unsigned long int seed, std::string simName) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.enable_varyenergy();
  in.set_simname(simName);
  in.set_spiralModel(SpiralModel::Steiman2010);

  runConfigured(in);
}

void runVarySlope(unsigned long int seed, std::string simName) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.enable_varyslope();
  in.set_simname(simName);
  in.set_spiralModel(SpiralModel::Steiman2010);

  runConfigured(in);
}

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    if (argc != 2) throw std::runtime_error("Usage: ./run params.ini");
    utils::Timer timer("timer for main");

    // runUniform(atoi(argv[1]), "test_solution");
    runFixed(atoi(argv[1]), 2. * cgs::kpc, "test_fixed2");
    runFixed(atoi(argv[1]), 4. * cgs::kpc, "test_fixed4");
    runFixed(atoi(argv[1]), 8. * cgs::kpc, "test_fixed8");
    runVaryEnergy(atoi(argv[1]), "test_varyesn");
    runVarySlope(atoi(argv[1]), "test_varyalpha");

  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
  }
  return EXIT_SUCCESS;
}