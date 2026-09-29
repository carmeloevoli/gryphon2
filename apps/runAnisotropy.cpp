#include <cstdlib>

#include "gryphon.h"

using namespace gryphon;

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "runAnisotropy");
    if (!args) return EXIT_SUCCESS;

    utils::Timer timer("timer for main");
    const auto input = utils::makeInput(*args);
    input.print();
    const core::RunOutput run(input, args->outdir);

    RandomNumberGenerator rng(input.seed());
    auto galaxy = galaxy::makeGalaxy(input);
    galaxy->generate(rng, false);
    const auto& events = galaxy->get_events();
    LOGD << "event size : " << events.size();

    auto kernel = kernel::makeGreenKernel(input);
    auto spectra = injection::makeInjectionSpectra(input, events, rng);
    core::CosmicRays cosmicRays(input, kernel, std::move(spectra), events);
    cosmicRays.run();
    core::dumpFluxAndDipole(run, cosmicRays);
  } catch (const std::exception& error) {
    LOGE << "exception caught with message: " << error.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
