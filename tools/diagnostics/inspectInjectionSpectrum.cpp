#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.set_injEmax(100. * cgs::TeV);
    in.set_efficiency(1.);
    in.set_injSlope(1.5);
    auto spectrum_15 = injection::SinglePowerLawSpectrum(in);
    in.set_injSlope(2.0);
    auto spectrum_20 = injection::SinglePowerLawSpectrum(in);
    in.set_injSlope(2.3);
    auto spectrum_23 = injection::SinglePowerLawSpectrum(in);
    auto energyAxis = utils::LogAxis<double>(1. * cgs::GeV, 1e6 * cgs::GeV, 1000);
    utils::OutputFile out("inspect_injection_spectrum.txt");
    out << "# E [GeV] - E2Q [erg]\n";
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / cgs::GeV << "\t";
      out << (pow2(E) * spectrum_15.get(E)) / cgs::erg << "\t";
      out << (pow2(E) * spectrum_20.get(E)) / cgs::erg << "\t";
      out << (pow2(E) * spectrum_23.get(E)) / cgs::erg << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}