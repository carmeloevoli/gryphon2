#include <iomanip>

#include "gryphon.h"

using namespace gryphon;

void dump_random_pulsars(const core::Input& in, size_t N) {
  RandomNumberGenerator rng(in.seed());

  utils::OutputFile out("inspect_random_pulsars.txt");
  out << "# P0 [ms] - B0 [G] - rotational energy [erg] - Emax [GeV] - tau0 [kyr]\n";
  out << std::scientific << std::setprecision(6);

  for (size_t i = 0; i < N; ++i) {
    const auto spectrum = injection::YoungPulsarsSpectrum(in, rng);
    out << spectrum.initialPeriod / cgs::msec << "\t";
    out << spectrum.surfaceMagneticField / cgs::gauss << "\t";
    out << spectrum.rotEnergy / cgs::erg << "\t";
    out << spectrum.Emax / cgs::GeV << "\t";
    out << spectrum.tau0 / cgs::kyr << "\t";
    out << "\n";
  }
}

void dump_injection_spectrum(const core::Input& in) {
  RandomNumberGenerator rng(in.seed());

  utils::OutputFile out("inspect_random_pulsars_spectrum.txt");
  out << "# E [GeV] - Q(E) [GeV^-1]\n";
  out << std::scientific << std::setprecision(6);

  const auto spectrum = injection::YoungPulsarsSpectrum(in, rng);

  LOGD << "YoungPulsarsSpectrum parameters:";
  LOGD << "  initial period: " << spectrum.initialPeriod / cgs::msec << " ms";
  LOGD << "  surface magnetic field: " << spectrum.surfaceMagneticField / cgs::gauss << " G";
  LOGD << "  maximum potential drop energy: " << spectrum.Emax / cgs::PeV << " PeV";
  LOGD << "  tau0: " << spectrum.tau0 / cgs::kyr << " kyr";
  LOGD << "  rotational energy: " << spectrum.rotEnergy / cgs::erg << " erg";
  LOGD << "  CR energy: " << spectrum.crEnergy / cgs::erg << " erg";

  const auto units = 1. / cgs::GeV;
  const auto energyAxis = utils::LogAxis<double>(cgs::TeV, 1e2 * cgs::PeV, 100);
  for (const auto& E : energyAxis) {
    out << E / cgs::GeV << "\t";
    out << spectrum.get(E) / units << "\t";
    out << "\n";
  }
}

int main() {
  try {
    utils::startup_information();
    auto in = core::Input();
    in.set_injectionModel(InjectionModel::YoungPulsars);
    in.set_efficiency(1.);
    in.print();

    constexpr size_t N = 100000;

    dump_random_pulsars(in, N);

    in.set_youngPulsarsRandomInitialPeriod(false);
    in.set_youngPulsarsRandomMagneticField(false);

    dump_injection_spectrum(in);

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
