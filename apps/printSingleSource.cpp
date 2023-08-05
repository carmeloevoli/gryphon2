#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.print();
    auto particle = std::make_shared<particle::FixedSpectrumParticle>(in);
    auto units = 1. / (1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr);
    const auto d = utils::Vector3d(cgs::kpc, 0., 0.);
    {
      auto timeAxis = utils::LogAxis<double>(cgs::kyr, 10. * cgs::Myr, 1000);
      utils::OutputFile out("test_single.txt");
      out << "# t [Myr] - flux []\n";
      out << std::scientific;
      for (auto t : timeAxis) {
        out << t / cgs::Myr << "\t";
        out << particle->get(cgs::TeV, t, d) / units << "\t";
        out << particle->get(10. * cgs::TeV, t, d) / units << "\t";
        out << particle->get(100. * cgs::TeV, t, d) / units << "\t";
        out << particle->get(cgs::PeV, t, d) / units << "\t";
        out << "\n";
      }
    }
    {
      auto energyAxis = utils::LogAxis<double>(cgs::TeV, cgs::PeV, 1000);
      auto timeAxis = utils::LogAxis<double>(0.01 * cgs::Myr, 0.1 * cgs::Myr, 10);
      utils::OutputFile out("test_single_spectrum.txt");
      out << "# E [GeV] - flux []\n";
      out << std::scientific;
      for (auto E : energyAxis) {
        out << E / cgs::GeV << "\t";
        for (auto t : timeAxis) out << particle->get(E, t, d) / units << "\t";
        out << "\n";
      }
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}