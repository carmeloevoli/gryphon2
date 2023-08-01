#include "gryphon.h"

using namespace gryphon;

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    if (argc != 2) throw std::runtime_error("Usage: ./run params.ini");
    utils::Timer timer("timer for main");

    // TODO check if output dir exists

    auto in = core ::Input(argv[1]);
    in.print();

    // if (rng->get_seed() == 0) input->save_on_file(rng);

    RandomNumberGenerator rng = utils::RNG<double>(in.seed);

    auto profile = core::SourceProfile(in);
    std::shared_ptr<galaxy::Galaxy> galaxy;
    switch (in.spiralModel) {
      case SpiralModel::Uniform:
        galaxy = std::make_shared<galaxy::GalaxyUniform>(in);
        break;
      default:
        throw std::invalid_argument("Spiral model not implemented yet");
    }
    galaxy->generate(rng);

    std::shared_ptr<particle::Particle> particle;
    switch (in.particleModel) {
      case ParticleModel::SingleSpectrum:
        particle = std::make_shared<particle::FixedSpectrumParticle>(in);
        break;
      default:
        throw std::invalid_argument("Particle model not implemented yet");
    }

    auto output = std::make_shared<core::OutputManager>(in);
    output->compute(galaxy, particle);
    output->dump();
  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
  }
  return EXIT_SUCCESS;
}