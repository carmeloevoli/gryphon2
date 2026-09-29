#include <cstdlib>
#include <stdexcept>

#include "gryphon.h"

using namespace gryphon;

namespace {

void runYoungPulsars(const core::Input& input, const core::RunOutput& run) {
  RandomNumberGenerator rng(input.seed());

  auto galaxyModel = galaxy::makeGalaxy(input);
  galaxyModel->generate(rng, false);
  const auto& events = galaxyModel->get_events();
  LOGD << "event size : " << events.size();

  auto kernel = kernel::makeGreenKernel(input);
  auto injectionSpectra = injection::makeInjectionSpectra(input, events, rng);
  core::CosmicRays cr(input, kernel, std::move(injectionSpectra), events);
  cr.run();

  core::dumpFluxAndDipole(run, cr);
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "runYoungPulsars");
    if (!args) return EXIT_SUCCESS;

    utils::Timer timer("timer for main");
    const auto input = utils::makeInput(*args);
    input.print();

    const core::RunOutput run(input, args->outdir);
    runYoungPulsars(input, run);

  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
