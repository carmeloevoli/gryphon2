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

void dumpOrigins(const core::Events& events,
                 const std::vector<galaxy::AssociationEventOrigin>& origins,
                 const galaxy::AssociationGenerationStats& stats,
                 const core::RunOutput& run) {
  if (events.size() != origins.size()) throw std::runtime_error("event/origin size mismatch");
  auto out = run.open("event_origins",
      "age [Myr] | x [kpc] | y [kpc] | z [kpc] | parent_id | birth_age [Myr] | "
      "delay [Myr] | centre_x [kpc] | centre_y [kpc] | centre_z [kpc]");
  out << std::scientific << std::setprecision(12);
  for (size_t i = 0; i < events.size(); ++i) {
    const auto& origin = origins[i];
    out << events[i]->age / cgs::Myr << "\t" << events[i]->pos / cgs::kpc << "\t"
        << origin.parentId << "\t" << origin.birthAge / cgs::Myr << "\t"
        << origin.delay / cgs::Myr << "\t" << origin.centre / cgs::kpc << "\n";
  }
  auto counts = run.open("population",
      "associations | field_explosions | clustered_explosions | retained_explosions");
  counts << stats.associations << "\t" << stats.fieldExplosions << "\t"
         << stats.clusteredExplosions << "\t" << stats.retainedExplosions << "\n";
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
    std::vector<galaxy::AssociationEventOrigin> origins;
    galaxyModel->generate(rng, false, input.syntheticAssociations() ? &origins : nullptr);

    const auto& events = galaxyModel->get_events();
    LOGD << "event size : " << events.size();

    const core::RunOutput run(input, args->outdir);
    dumpEvents(events, run);
    if (input.syntheticAssociations())
      dumpOrigins(events, origins, galaxyModel->association_stats(), run);

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
