#include <iomanip>

#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
    auto in = core::Input();
    in.set_injectionModel(InjectionModel::RandomEmax);
    in.set_efficiency(1.);
    in.set_injSlope(2.32);
    in.set_injEmax(5. * cgs::PeV);
    in.set_injEmaxSigmaDex(0.4);
    in.set_injEmaxMin(100. * cgs::TeV);
    in.set_injEmaxMax(100. * cgs::PeV);
    in.print();

    constexpr size_t N = 10000;
    RandomNumberGenerator rng(in.seed());

    utils::OutputFile out("inspect_random_emax.txt");
    out << "# Emax [TeV]\n";
    out << std::scientific << std::setprecision(6);
    for (size_t i = 0; i < N; ++i) {
      const auto spectrum = injection::RandomEmaxSpectrum(in, rng);
      out << spectrum.crEnergy << "\t";
      out << spectrum.Emax / cgs::TeV << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
