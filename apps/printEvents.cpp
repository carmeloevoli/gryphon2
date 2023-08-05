#include "gryphon.h"

using namespace gryphon;

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
      auto events = galaxy->get_events();

      utils::OutputFile out("test_profile.txt");
      out << "#\n";
      out << std::scientific << std::setprecision(6);
      for (auto& event : events) {
        out << event->age / cgs::Myr << "\t";
        out << event->pos / cgs::kpc << "\t";
        out << event->pos.getModule() / cgs::kpc << "\t";
        out << "\n";
      }
    }
    {
      auto galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
      galaxy->generate(rng, false);
      auto events = galaxy->get_events();

      utils::OutputFile out("test_pattern.txt");
      out << "#\n";
      out << std::scientific << std::setprecision(6);
      for (auto& event : events) {
        out << event->age / cgs::Myr << "\t";
        out << event->pos / cgs::kpc << "\t";
        out << event->pos.getModule() / cgs::kpc << "\t";
        out << "\n";
      }
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}