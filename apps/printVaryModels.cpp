#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.set_efficiency(1.);
    in.set_injEmax(1e3 * cgs::TeV);
    in.enable_varyenergy();
    in.enable_varyslope();
    in.print();

    size_t N = 10000;
    RandomNumberGenerator rng = utils::RNG<double>(in.seed);

    utils::OutputFile out("test_random_spectrum.txt");
    out << "# CR energy [erg] - slope \n";
    out << std::scientific << std::setprecision(6);
    for (size_t i = 0; i < N; ++i) {
      auto spectrum = injection::GalacticRandomSpectrum(in, rng);
      out << spectrum.crEnergy / cgs::erg << "\t";
      out << spectrum.alpha << "\t";
      out << spectrum.Q0 / spectrum.crEnergy << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}