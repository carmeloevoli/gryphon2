// Samples the young-pulsar birth model of a configuration and dumps both the
// resulting distribution of birth properties and the injection spectrum of the
// reference pulsar (randomness switched off).
#include <cstdlib>
#include <iomanip>

#include "gryphon.h"

using namespace gryphon;

namespace {

constexpr size_t kSamples = 100000;

void dumpBirthProperties(const core::Input& in, const core::RunOutput& run) {
  RandomNumberGenerator rng(in.seed());

  auto out = run.open("pulsars",
                      "P0 [ms] | B0 [G] | rotational energy [erg] | Emax [GeV] | tau0 [kyr]");
  out << std::scientific << std::setprecision(6);

  for (size_t i = 0; i < kSamples; ++i) {
    const auto spectrum = injection::YoungPulsarsSpectrum(in, rng);
    out << spectrum.initialPeriod / cgs::msec << "\t";
    out << spectrum.surfaceMagneticField / cgs::gauss << "\t";
    out << spectrum.rotEnergy / cgs::erg << "\t";
    out << spectrum.Emax / cgs::GeV << "\t";
    out << spectrum.tau0 / cgs::kyr << "\n";
  }
}

void dumpReferenceSpectrum(const core::Input& in, const core::RunOutput& run) {
  RandomNumberGenerator rng(in.seed());
  const auto spectrum = injection::YoungPulsarsSpectrum(in, rng);

  LOGD << "reference YoungPulsarsSpectrum:";
  LOGD << "  initial period: " << spectrum.initialPeriod / cgs::msec << " ms";
  LOGD << "  surface magnetic field: " << spectrum.surfaceMagneticField / cgs::gauss << " G";
  LOGD << "  maximum potential drop energy: " << spectrum.Emax / cgs::PeV << " PeV";
  LOGD << "  tau0: " << spectrum.tau0 / cgs::kyr << " kyr";
  LOGD << "  rotational energy: " << spectrum.rotEnergy / cgs::erg << " erg";
  LOGD << "  CR energy: " << spectrum.crEnergy / cgs::erg << " erg";

  auto out = run.open("injection", "E [GeV] | Q(E) [GeV^-1]");
  out << std::scientific << std::setprecision(6);

  const auto units = 1. / cgs::GeV;
  const auto energyAxis = utils::LogAxis<double>(cgs::TeV, 1e2 * cgs::PeV, 100);
  for (const auto& E : energyAxis) {
    out << E / cgs::GeV << "\t";
    out << spectrum.get(E) / units << "\n";
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "inspectRandomPulsars");
    if (!args) return EXIT_SUCCESS;

    auto input = utils::makeInput(*args);
    input.set_injectionModel(InjectionModel::YoungPulsars);
    input.print();

    const core::RunOutput run(input, args->outdir);
    dumpBirthProperties(input, run);

    // The reference spectrum is the one of a pulsar sitting at the centre of
    // the birth distributions, so the draws are switched off.
    input.set_youngPulsarsRandomInitialPeriod(false);
    input.set_youngPulsarsRandomMagneticField(false);
    dumpReferenceSpectrum(input, run);

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
