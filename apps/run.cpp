#include "gryphon.h"

using namespace gryphon;

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    if (argc != 2) throw std::runtime_error("Usage: ./run params.ini");
    utils::Timer timer("timer for main");

    // TODO check if output dir exists

    auto in = core ::Input();
    in.set_seed(atoi(argv[1]));
    in.set_maxtime(0.1 * cgs::Myr);
    in.set_simname("testNew");
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

    particle::Particles particles;
    particles.reserve(galaxy->size());

    auto events = galaxy->get_events();
    for (auto& event : events) {
      auto particle = std::make_shared<particle::FixedSpectrumParticle>(in, event, rng);
      particles.emplace_back(particle);
    }

    LOGD << "particle size : " << particles.size();

    // std::shared_ptr<particle::Particle> particle;
    // switch (in.particleModel) {
    //   case ParticleModel::FixedSpectrum:
    //     particle = std::make_shared<particle::FixedSpectrumParticle>(in);
    //     break;
    //   default:
    //     throw std::invalid_argument("Particle model not implemented yet");
    // }

    auto output = std::make_shared<core::OutputManager>(in);
    output->compute(particles);
    output->dump();
  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
  }
  return EXIT_SUCCESS;
}