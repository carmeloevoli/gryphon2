// Dumps the source events of one realization: ages and positions, as generated
// by the spiral model of the configuration.
#include <cstdlib>
#include <iomanip>

#include "gryphon.h"

using namespace gryphon;

namespace {

void dumpEvents(const core::Events& events, const core::RunOutput& run) {
  const double invMyr = 1. / cgs::Myr;
  const double invKpc = 1. / cgs::kpc;

  auto out = run.open(
      "events", "age [Myr] | x [kpc] | y [kpc] | z [kpc] | distance from Sun [kpc]");
  out << std::scientific << std::setprecision(6);

  for (const auto& event : events) {
    if (!event) continue;
    out << event->age * invMyr << "\t";
    out << event->pos * invKpc << "\t";
    out << event->pos.getModule() * invKpc << "\n";
  }
}


}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "inspectEvents");
    if (!args) return EXIT_SUCCESS;

    const auto input = utils::makeInput(*args);
    input.print();

    RandomNumberGenerator rng(input.seed());
    auto galaxyModel = galaxy::makeGalaxy(input);
    galaxyModel->generate(rng, false);

    const auto& events = galaxyModel->get_events();
    LOGD << "event size : " << events.size();

    const core::RunOutput run(input, args->outdir);
    dumpEvents(events, run);

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
