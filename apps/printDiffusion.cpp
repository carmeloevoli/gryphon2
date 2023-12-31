#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.set_halosize(4. * cgs::kpc);
    in.print();
    auto D = core::DiffusionCoefficient(in);
    auto energyAxis = utils::LogAxis<double>(1e2 * cgs::GeV, 1e6 * cgs::GeV, 100);
    utils::OutputFile out("test_diffusion.txt");
    out << "# E [GeV] - D [cm2/s]\n";
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / cgs::GeV << "\t";
      out << D.get(E) / (cgs::kpc2 / cgs::Myr) << "\t";
      out << D.escapeTimescale(E) / cgs::Myr << "\t";
      out << "\n";
    }
    LOGI << "tau : " << D.escapeTimescale(cgs::TeV) / cgs::Myr << " Myr\n";
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}