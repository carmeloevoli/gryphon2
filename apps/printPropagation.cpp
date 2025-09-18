#include "gryphon.h"

using namespace gryphon;

void printTimescales(const core::Propagation& prop) {
  auto energyAxis = utils::LogAxis<double>(cgs::GeV, 1e4 * cgs::GeV, 100);
  utils::OutputFile out("test_timescales.txt");
  out << "# E [GeV] - D [cm2/s]\n";
  out << std::scientific;
  for (auto E : energyAxis) {
    out << E / cgs::GeV << "\t";
    out << pow2(prop.H) / 2. / prop.D(E) / cgs::Myr << "\t";
    out << E / prop.b(E) / cgs::Myr << "\t";
    out << "\n";
  }
}

void printEstar(const core::Propagation& prop) {
  auto units = cgs::GeV;
  auto energyAxis = utils::LogAxis<double>(cgs::GeV, 1e4 * cgs::GeV, 100);
  utils::OutputFile out("test_msp_Estar.txt");
  out << "# E [GeV] - E_star []\n";
  out << std::scientific;
  for (auto E : energyAxis) {
    out << E / cgs::GeV << "\t";
    out << prop.Estar(E, 1e1 * cgs::kyr) / units << "\t";
    out << prop.Estar(E, 1e2 * cgs::kyr) / units << "\t";
    out << prop.Estar(E, 1e3 * cgs::kyr) / units << "\t";
    out << prop.Estar(E, 1e4 * cgs::kyr) / units << "\t";
    out << prop.Estar(E, 1e5 * cgs::kyr) / units << "\t";
    out << "\n";
  }
}

void printLambda(const core::Propagation& prop) {
  auto energyAxis = utils::LogAxis<double>(cgs::GeV, cgs::PeV, 10000);
  utils::OutputFile out("test_msp_lambda.txt");
  out << "# E [GeV] - E_star []\n";
  out << std::scientific;
  for (auto E : energyAxis) {
    out << E / cgs::GeV << "\t";
    out << std::sqrt(prop.lambda2(E, 1e10 * E)) / cgs::kpc << "\t";
    out << prop.tau(E, 1e10 * E) / cgs::Myr << "\t";
    out << "\n";
  }
}

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.set_refEnergy(cgs::GeV);
    in.set_D0_over_H(0.35e28 * cgs::cm2 / cgs::sec / cgs::kpc);
    in.set_delta(0.56);
    in.set_halosize(5. * cgs::kpc);
    in.set_galaxyRadius(20. * cgs::kpc);
    in.set_Bfield(5. * cgs::microgauss);
    in.print();

    auto prop = core::Propagation(in);
    printTimescales(prop);
    printEstar(prop);
    printLambda(prop);
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}