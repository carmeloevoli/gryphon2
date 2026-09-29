// Per-realization participation number, using the production galaxy,
// injection factory, and Green function without changing production outputs.
#include <cstdint>
#include <cstdlib>
#include <iomanip>

#include "gryphon.h"
#include "gryphon/core/participation.h"

using namespace gryphon;

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "inspectEffectiveSources");
    if (!args) return EXIT_SUCCESS;
    const auto input = utils::makeInput(*args);
    input.print();

    RandomNumberGenerator rng(input.seed());
    auto galaxy = galaxy::makeGalaxy(input);
    galaxy->generate(rng, false);
    const auto& events = galaxy->get_events();

    // FNV-1a fingerprint of the ordered ages/positions in native double bytes.
    // Only an identity check between these matched runs on the same platform,
    // not a portable catalogue format or cryptographic content digest.
    std::uint64_t fingerprint = 14695981039346656037ULL;
    size_t count = 0;
    for (const auto& event : events) {
      if (!event) continue;
      ++count;
      for (const double value : {event->age, event->pos.x, event->pos.y, event->pos.z}) {
        const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
        for (size_t j = 0; j < sizeof(value); ++j) {
          fingerprint ^= bytes[j];
          fingerprint *= 1099511628211ULL;
        }
      }
    }
    // Same ordering and RNG state as runAnisotropy: first catalogue, then
    // random injection properties. The diagnostic introduces no random draws.
    const auto spectra = injection::makeInjectionSpectra(input, events, rng);
    const auto kernel = kernel::makeGreenKernel(input);
    const auto energy = utils::LogAxis<double>(input.E_min(), input.E_max(), input.E_size());
    const core::RunOutput run(input, args->outdir);
    auto out = run.open("participation",
        "E [GeV] | I [GeV^-1 m^-2 s^-1 sr^-1] | N_eff | q_max");
    out << "# source_count = " << count << "\n";
    out << "# catalogue_fnv1a64 = " << std::hex << fingerprint << std::dec << "\n";
    out << std::scientific << std::setprecision(16);
    const double flux_units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;
    for (const double E : energy) {
      double flux = 0.;
      core::Participation weights;
      for (size_t j = 0; j < events.size(); ++j) {
        if (!events[j] || !spectra[j]) continue;
        const auto* spectrum = spectra[j].get();
        const double contribution = kernel->contribution(
            E, events[j]->age, events[j]->pos,
            [spectrum](double Eprime) { return spectrum->get(Eprime); }).flux;
        flux += contribution;
        weights.add(contribution);
      }
      out << E / cgs::GeV << '\t' << flux / flux_units << '\t'
          << weights.effectiveCount() << '\t' << weights.maxFraction() << '\n';
    }
  } catch (const std::exception& error) {
    LOGE << "exception caught with message: " << error.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
