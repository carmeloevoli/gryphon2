#include "gryphon.h"

using namespace gryphon;

namespace {

void runAndDump(const core::Input& in, const core::Events& events, RandomNumberGenerator& rng) {
  auto kernel = kernel::makeGreenKernel(in);
  auto injectionSpectrum = injection::makeInjectionSpectrum(in, rng);
  core::CosmicRays cr(in, kernel, injectionSpectrum, events);
  cr.run();

  const double units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;
  utils::OutputFile out(in.simname + "_" + std::to_string(in.seed) + ".txt");
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

void runUniform(unsigned long int seed, std::string simname) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.set_simname(simname);
  // in.print();

  RandomNumberGenerator rng = utils::RNG<double>(in.seed);

  auto galaxy = std::make_shared<galaxy::GalaxyUniform>(in);
  galaxy->generate(rng, false);

  auto events = galaxy->get_events();
  runAndDump(in, events, rng);
}

void runFixed(unsigned long int seed, double haloSize, std::string simName) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.set_halosize(haloSize);
  in.set_simname(simName);
  // in.print();

  RandomNumberGenerator rng = utils::RNG<double>(in.seed);

  auto galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
  galaxy->generate(rng, false);

  auto events = galaxy->get_events();
  runAndDump(in, events, rng);
}

void runVaryEnergy(unsigned long int seed, std::string simName) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.enable_varyenergy();
  in.set_simname(simName);

  RandomNumberGenerator rng = utils::RNG<double>(in.seed);

  auto galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
  galaxy->generate(rng, false);

  auto events = galaxy->get_events();
  runAndDump(in, events, rng);
}

void runVarySlope(unsigned long int seed, std::string simName) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.enable_varyslope();
  in.set_simname(simName);

  RandomNumberGenerator rng = utils::RNG<double>(in.seed);

  auto galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
  galaxy->generate(rng, false);

  auto events = galaxy->get_events();
  runAndDump(in, events, rng);
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