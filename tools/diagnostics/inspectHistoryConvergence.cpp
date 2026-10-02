// Compare nested age cuts of exactly the same catalogue and injection draws.
#include <cstdlib>
#include <iomanip>
#include "gryphon.h"

using namespace gryphon;

int main(int argc, char* argv[]) {
  try {
    const auto args = utils::parseCommandLine(argc, argv, "inspectHistoryConvergence");
    if (!args) return EXIT_SUCCESS;
    const auto input = utils::makeInput(*args);
    if (input.max_time() < 300. * cgs::Myr)
      throw std::invalid_argument("history diagnostic requires at least 300 Myr");
    RandomNumberGenerator rng(input.seed());
    auto galaxy = galaxy::makeGalaxy(input);
    galaxy->generate(rng, false);
    const auto& events = galaxy->get_events();
    const auto spectra = injection::makeInjectionSpectra(input, events, rng);
    const auto kernel = kernel::makeGreenKernel(input);
    const core::RunOutput run(input, args->outdir);
    auto out = run.open("history", "E [GeV] | I_100 | I_200 | I_300");
    out << std::scientific << std::setprecision(16);
    const double units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;
    for (double E : utils::LogAxis<double>(input.E_min(), input.E_max(), input.E_size())) {
      double young = 0., middle = 0., old = 0.;
      for (size_t j = 0; j < events.size(); ++j) {
        const auto* spectrum = spectra[j].get();
        const double value = kernel->contribution(E, events[j]->age, events[j]->pos,
            [spectrum](double e) { return spectrum->get(e); }).flux;
        if (events[j]->age <= 100. * cgs::Myr) young += value;
        else if (events[j]->age <= 200. * cgs::Myr) middle += value;
        else if (events[j]->age <= 300. * cgs::Myr) old += value;
      }
      out << E / cgs::GeV << '\t' << young / units << '\t'
          << (young + middle) / units << '\t' << (young + middle + old) / units << '\n';
    }
  } catch (const std::exception& error) {
    LOGE << error.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
