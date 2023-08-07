#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.set_efficiency(1.);
    in.enable_varyenergy();
    in.enable_varyslope();
    in.print();

    RandomNumberGenerator rng = utils::RNG<double>(in.seed);
    auto event = std::make_shared<core::Event>();
    size_t N = 10000;

    utils::OutputFile out("test_vary.txt");
    out << "# CR energy [erg] - slope \n";
    out << std::scientific << std::setprecision(6);
    for (size_t i = 0; i < N; ++i) {
      auto particle = particle::FixedSpectrumParticle(in, event, rng);
      out << particle.crEnergy / cgs::erg << "\t";
      out << particle.alpha << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}