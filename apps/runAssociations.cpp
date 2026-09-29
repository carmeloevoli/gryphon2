#include <cstdlib>

#include "gryphon.h"

using namespace gryphon;

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "runAssociations");
    if (!args) return EXIT_SUCCESS;
    const auto input = utils::makeInput(*args);
    if (!input.syntheticAssociations())
      throw std::invalid_argument("runAssociations requires syntheticassociations = true");
    const core::RunOutput run(input, args->outdir);
    RandomNumberGenerator rng(input.seed());
    auto population = galaxy::makeGalaxy(input);
    population->generate(rng, false);
    const auto& events = population->get_events();
    auto kernel = kernel::makeGreenKernel(input);
    auto spectra = injection::makeInjectionSpectra(input, events, rng);
    core::CosmicRays cosmicRays(input, kernel, std::move(spectra), events);
    cosmicRays.run();
    core::dumpFluxAndDipole(run, cosmicRays);
    const auto& stats = population->association_stats();
    auto out = run.open("population", "parents | field_SN | clustered_SN | retained_SN | expected_SN");
    out << stats.associations << '\t' << stats.fieldExplosions << '\t'
        << stats.clusteredExplosions << '\t' << stats.retainedExplosions << '\t'
        << input.sn_rate() * input.max_time() << '\n';
  } catch (const std::exception& error) {
    LOGE << error.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
