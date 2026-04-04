#include "gryphon.h"

using namespace gryphon;

namespace {

void dumpEvents(const core::Events& events, const std::string& filename) {
  const double invMyr = 1. / cgs::Myr;
  const double invKpc = 1. / cgs::kpc;

  utils::OutputFile out(filename);
  out << "#\n";
  out << std::scientific << std::setprecision(6);
  for (const auto& event : events) {
    out << event->age * invMyr << "\t";
    out << event->pos * invKpc << "\t";
    out << event->pos.getModule() * invKpc << "\t";
    out << "\n";
  }
}

}  // namespace

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.set_maxtime(cgs::Myr);
    in.print();

    RandomNumberGenerator rng = utils::RNG<double>(in.seed);
    {
      auto galaxy = std::make_shared<galaxy::GalaxyJelly>(in);
      galaxy->generate(rng, false);
      const auto& events = galaxy->get_events();
      dumpEvents(events, "test_galaxy_jelly.txt");
    }
    {
      auto galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
      galaxy->generate(rng, false);
      const auto& events = galaxy->get_events();
      dumpEvents(events, "test_galaxy_spirals.txt");
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}