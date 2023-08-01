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
    auto D = core::DiffusionCoefficient(in);
    auto profile = core::SourceProfile(in);

    std::shared_ptr<galaxy::Galaxy> galaxy;
    switch (in.spiralModel) {
      case SpiralModel::Uniform:
        galaxy = std::make_shared<galaxy::GalaxyUniform>(in);
        break;
      default:
        throw std::runtime_error("Spiral model not implemented yet");
    }
    galaxy->generate(rng);

    // std::shared_ptr<GRYPHON::Galaxy> galaxy;
    // switch (input->spiralModel) {
    //   case GRYPHON::SpiralModel::Uniform:
    //     galaxy = std::make_shared<GRYPHON::GalaxyUniform>(input, rng);
    //     break;
    //   case GRYPHON::SpiralModel::Jelly:
    //     galaxy = std::make_shared<GRYPHON::GalaxyJelly>(input, rng, profile);
    //     break;
    //   case GRYPHON::SpiralModel::Steiman2010:
    //     galaxy = std::make_shared<GRYPHON::GalaxySteiman2010>(input, rng, profile);
    //     break;
    //   default:
    //     throw std::runtime_error("Spiral model not implemented yet");
    // }

    // std::shared_ptr<GRYPHON::Particle> particle;
    // switch (input->particleType) {
    //   case GRYPHON::ParticleType::Protons:
    //     particle = std::make_shared<GRYPHON::Nuclei>(input, Db, GRYPHON::PID(1, 1));
    //     break;
    //   case GRYPHON::ParticleType::Helium:
    //     particle = std::make_shared<GRYPHON::Nuclei>(input, Db, GRYPHON::PID(2, 4));
    //     break;
    //   case GRYPHON::ParticleType::Electrons:
    //     particle = std::make_shared<GRYPHON::Electrons>(input, Db);
    //     break;
    //   case GRYPHON::ParticleType::Positrons:
    //     particle = std::make_shared<GRYPHON::Positrons>(input, Db);
    //     break;
    //   case GRYPHON::ParticleType::PositronsRunning:
    //     particle = std::make_shared<GRYPHON::PositronsRunning>(input, Db);
    //     break;
    //   default:
    //     throw std::runtime_error("Particle model not implemented yet");
    // }

    // auto output = std::make_shared<GRYPHON::OutputManager>(particle, input);

    // if (input->runMode == GRYPHON::RunMode::Flux) {
    //   galaxy->generate_events();
    //   particle->run_flux(galaxy);
    //   output->save_flux();
    // }

    // if (input->runMode == GRYPHON::RunMode::Anisotropy) {
    //   galaxy->generate_events();
    //   particle->run_anisotropy(galaxy);
    //   output->save_anisotropy();
    // }
  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
  }
  return EXIT_SUCCESS;
}