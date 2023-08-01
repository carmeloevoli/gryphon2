#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.print();
    auto particle = std::make_shared<particle::FixedSpectrumParticle>(in);
    auto energyAxis = utils::LogAxis<double>(1e1 * cgs::GeV, 1e5 * cgs::GeV, 400);
    auto age = 0.5 * cgs::Myr;
    auto units = 1. / (1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr);
    utils::OutputFile out("test_single.txt");
    out << "# E [GeV] - flux []\n";
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / cgs::GeV << "\t";
      out << pow2(E) * particle->Q(E) / cgs::erg << "\t";
      out << particle->get(E, age, utils::Vector3d(0.3 * cgs::kpc, 0., 0.)) / units << "\t";
      out << particle->get(E, age, utils::Vector3d(0.5 * cgs::kpc, 0., 0.)) / units << "\t";
      out << particle->get(E, age, utils::Vector3d(1. * cgs::kpc, 0., 0.)) / units << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}