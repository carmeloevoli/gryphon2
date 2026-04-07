#include <iomanip>

#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
    auto in = core::Input();
    in.set_injectionModel(InjectionModel::PWN);
    in.set_efficiency(1.);
    in.set_pwnP0(100. * cgs::msec);
    in.set_pwnSigmaP0(30. * cgs::msec);
    in.set_pwnAlpha1(1.5);
    in.set_pwnAlpha2(2.6);
    in.set_pwnEmin(1. * cgs::GeV);
    in.set_pwnEbreak(100. * cgs::GeV);
    in.print();

    constexpr size_t N = 100000;
    RandomNumberGenerator rng(in.seed());

    utils::OutputFile out("inspect_random_pulsars.txt");
    out << "# P0 [ms] - vx [km/sec] - vy [km/sec] - vz [km/sec] - |v| [km/sec]\n";
    out << std::scientific << std::setprecision(6);
    for (size_t i = 0; i < N; ++i) {
      const auto spectrum = injection::PWNSpectrum(in, rng);
      out << spectrum.initialPeriod / cgs::msec << "\t";
      out << spectrum.kickVelocity.getX() / (cgs::km / cgs::sec) << "\t";
      out << spectrum.kickVelocity.getY() / (cgs::km / cgs::sec) << "\t";
      out << spectrum.kickVelocity.getZ() / (cgs::km / cgs::sec) << "\t";
      out << spectrum.kickSpeed / (cgs::km / cgs::sec) << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
