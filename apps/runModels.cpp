#include "gryphon.h"

using namespace gryphon;

void runUniform(unsigned long int seed, std::string simname) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.set_simname(simname);
  // in.print();

  RandomNumberGenerator rng = utils::RNG<double>(in.seed);

  auto galaxy = std::make_shared<galaxy::GalaxyUniform>(in);
  galaxy->generate(rng);

  particle::Particles particles;
  particles.reserve(galaxy->size());

  auto events = galaxy->get_events();
  for (auto& event : events) {
    auto particle = std::make_shared<particle::FixedSpectrumParticle>(in, event, rng);
    particles.emplace_back(particle);
  }

  auto output = std::make_shared<core::OutputManager>(in);
  output->compute(particles);
  output->dump();
}

void runFixed(unsigned long int seed, double haloSize, std::string simName) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.set_halosize(haloSize);
  in.set_simname(simName);
  // in.print();

  RandomNumberGenerator rng = utils::RNG<double>(in.seed);

  auto galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
  galaxy->generate(rng);

  particle::Particles particles;
  particles.reserve(galaxy->size());

  auto events = galaxy->get_events();
  for (auto& event : events) {
    auto particle = std::make_shared<particle::FixedSpectrumParticle>(in, event, rng);
    particles.emplace_back(particle);
  }

  auto output = std::make_shared<core::OutputManager>(in);
  output->compute(particles);
  output->dump();
}

void runVaryEnergy(unsigned long int seed, std::string simName) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.enable_varyenergy();
  in.set_simname(simName);

  RandomNumberGenerator rng = utils::RNG<double>(in.seed);

  auto galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
  galaxy->generate(rng);

  particle::Particles particles;
  particles.reserve(galaxy->size());

  auto events = galaxy->get_events();
  for (auto& event : events) {
    auto particle = std::make_shared<particle::FixedSpectrumParticle>(in, event, rng);
    particles.emplace_back(particle);
  }

  auto output = std::make_shared<core::OutputManager>(in);
  output->compute(particles);
  output->dump();
}

void runVarySlope(unsigned long int seed, std::string simName) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.enable_varyslope();
  in.set_simname(simName);

  RandomNumberGenerator rng = utils::RNG<double>(in.seed);

  auto galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
  galaxy->generate(rng);

  particle::Particles particles;
  particles.reserve(galaxy->size());

  auto events = galaxy->get_events();
  for (auto& event : events) {
    auto particle = std::make_shared<particle::FixedSpectrumParticle>(in, event, rng);
    particles.emplace_back(particle);
  }

  auto output = std::make_shared<core::OutputManager>(in);
  output->compute(particles);
  output->dump();
}

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    if (argc != 2) throw std::runtime_error("Usage: ./run params.ini");
    utils::Timer timer("timer for main");

    // runUniform(atoi(argv[1]), "test_solution");
    // runFixed(atoi(argv[1]), 2. * cgs::kpc, "test_fixed2");
    // runFixed(atoi(argv[1]), 4. * cgs::kpc, "test_fixed4");
    // runFixed(atoi(argv[1]), 8. * cgs::kpc, "test_fixed8");
    // runVaryEnergy(atoi(argv[1]), "test_varyesn");
    runVarySlope(atoi(argv[1]), "test_varyalpha");

  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
  }
  return EXIT_SUCCESS;
}