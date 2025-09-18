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
    in.set_simEmin(0.1 * cgs::TeV);
    in.set_simEmax(1e3 * cgs::TeV);
    in.set_simEsize(4 * 16);
    in.set_maxtime(100 * cgs::Myr);
    in.set_halosize(2 * cgs::kpc);
    in.set_simname("base");
    in.print();

    utils::OutputFile out("test_NSources.txt");
    out << "# E [GeV] - N\n";
    out << std::scientific;

    // Compute the number of sources for each energy bin
    auto EAxis = utils::LogAxis<double>(in.E_min, in.E_max, 5);
    for (auto E : EAxis) {
      // Compute the number of sources for the given energy
      double D_H = in.D0_over_H * std::pow(E / in.E_0, in.delta);
      double N = 0.5 * in.sn_rate * pow3(in.H) / pow2(in.R_g) / D_H;
      out << E / cgs::GeV << "\t" << int(N) << "\n";
    }
    // for (auto z : zAxis) {
    //   out << z / cgs::kpc << "\t";
    //   out << utils::halo_function(pow2(0.1 * in.H), in.H, z, 0.) << "\t";
    //   out << utils::halo_function(pow2(0.5 * in.H), in.H, z, 0.) << "\t";
    //   out << utils::halo_function(pow2(in.H), in.H, z, 0.) << "\t";
    //   out << utils::halo_function(pow2(3. * in.H), in.H, z, 0.) << "\t";
    //   out << "\n";
    //}

    // RandomNumberGenerator rng = utils::RNG<double>(in.seed);

    // std::shared_ptr<galaxy::Galaxy> galaxy;
    // switch (in.spiralModel) {
    //   case SpiralModel::Uniform:
    //     galaxy = std::make_shared<galaxy::GalaxyUniform>(in);
    //     break;
    //   case SpiralModel::Jelly:
    //     galaxy = std::make_shared<galaxy::GalaxyJelly>(in);
    //     break;
    //   case SpiralModel::Steiman2010:
    //     galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
    //     break;
    //   default:
    //     throw std::invalid_argument("Spiral model not implemented yet");
    // }
    // galaxy->generate(rng);

    // particle::Particles particles;
    // particles.reserve(galaxy->size());

    // auto events = galaxy->get_events();
    // for (auto& event : events) {
    //   auto particle = std::make_shared<particle::FixedSpectrumParticle>(in, event, rng);
    //   particles.emplace_back(particle);
    // }

    // LOGD << "particle size : " << particles.size();

    // // std::shared_ptr<particle::Particle> particle;
    // // switch (in.particleModel) {
    // //   case ParticleModel::FixedSpectrum:
    // //     particle = std::make_shared<particle::FixedSpectrumParticle>(in);
    // //     break;
    // //   default:
    // //     throw std::invalid_argument("Particle model not implemented yet");
    // // }

    // auto output = std::make_shared<core::OutputManager>(in);
    // output->compute(particles);
    // output->dump();
  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
  }
  return EXIT_SUCCESS;
}