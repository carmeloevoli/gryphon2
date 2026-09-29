// Samples the source-to-source scatter of the GalacticRandom injection model.
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
    const auto args = utils::parseCommandLine(argc, argv, "inspectRandomInjection");
    if (!args) return EXIT_SUCCESS;

    auto input = utils::makeInput(*args);
    input.set_injectionModel(InjectionModel::GalacticRandom);
    input.print();

    RandomNumberGenerator rng(input.seed());

    const core::RunOutput run(input, args->outdir);
    auto out = run.open("randominjection", "crEnergy [erg] | slope | Q0/crEnergy [GeV^-1]");
    out << std::scientific << std::setprecision(6);

    for (size_t i = 0; i < kSamples; ++i) {
      const auto spectrum = injection::GalacticRandomSpectrum(input, rng);
      out << spectrum.crEnergy / cgs::erg << "\t";
      out << spectrum.alpha << "\t";
      out << spectrum.Q0 / spectrum.crEnergy << "\n";
    }

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
