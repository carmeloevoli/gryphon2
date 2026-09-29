// Samples the log-normal cutoff-energy model of a configuration.
#include <cstdlib>
#include <iomanip>

#include "gryphon.h"

using namespace gryphon;

namespace {

constexpr size_t kSamples = 10000;

}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "inspectRandomEmax");
    if (!args) return EXIT_SUCCESS;

    auto input = utils::makeInput(*args);
    input.set_injectionModel(InjectionModel::RandomEmax);
    input.print();

    RandomNumberGenerator rng(input.seed());

    const core::RunOutput run(input, args->outdir);
    auto out = run.open("randomemax", "crEnergy [erg] | Emax [TeV]");
    out << std::scientific << std::setprecision(6);

    for (size_t i = 0; i < kSamples; ++i) {
      const auto spectrum = injection::RandomEmaxSpectrum(input, rng);
      out << spectrum.crEnergy / cgs::erg << "\t";
      out << spectrum.Emax / cgs::TeV << "\n";
    }

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
