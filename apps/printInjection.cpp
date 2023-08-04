#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.set_Emax(100. * cgs::TeV);
    in.set_efficiency(1.);
    in.print();
    auto particle = particle::FixedSpectrumParticle(in);
    auto energyAxis = utils::LogAxis<double>(1e2 * cgs::GeV, 1e6 * cgs::GeV, 100);
    utils::OutputFile out("test_injection.txt");
    out << "# E [GeV] - D [cm2/s]\n";
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / cgs::GeV << "\t";
      out << (pow2(E) * particle.Q(E)) / cgs::erg << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}