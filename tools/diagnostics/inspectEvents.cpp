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

void generateAndDump(core::Input& in, SpiralModel model, RandomNumberGenerator& rng,
                     const std::string& filename) {
  in.set_spiralModel(model);
  auto galaxyModel = galaxy::makeGalaxy(in);
  galaxyModel->generate(rng, false);
  dumpEvents(galaxyModel->get_events(), filename);
}

}  // namespace

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.set_maxtime(cgs::Myr);
    in.set_rate(1. / 59. / cgs::year);
    in.print();

    RandomNumberGenerator rng(in.seed());
    generateAndDump(in, SpiralModel::Jelly, rng, "inspect_events_jelly.txt");
    generateAndDump(in, SpiralModel::Steiman2010, rng, "inspect_events_steiman2010.txt");
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}
