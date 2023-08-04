#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.set_maxtime(1. * cgs::Myr);
    in.print();

    RandomNumberGenerator rng = utils::RNG<double>(in.seed);

    std::shared_ptr<galaxy::Galaxy> galaxy;
    switch (in.spiralModel) {
      case SpiralModel::Uniform:
        galaxy = std::make_shared<galaxy::GalaxyUniform>(in);
        break;
      case SpiralModel::Jelly:
        galaxy = std::make_shared<galaxy::GalaxyJelly>(in);
        break;
      case SpiralModel::Steiman2010:
        galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
        break;
      default:
        throw std::runtime_error("Spiral model not implemented yet");
    }
    galaxy->generate(rng, false);
    auto events = galaxy->get_events();

    utils::OutputFile out("test_events.txt");
    out << "#\n";
    out << std::scientific << std::setprecision(6);
    for (auto& event : events) {
      out << event->age / cgs::Myr << "\t";
      out << event->pos / cgs::kpc << "\t";
      out << event->pos.getModule() / cgs::kpc << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}