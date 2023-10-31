#include "gryphon.h"

using namespace gryphon;

void printLis() {
  auto H_lis = Boschini2017::H_LIS();
  auto He_lis = Boschini2017::He_LIS();
  auto energyAxis = utils::LogAxis<double>(0.1 * cgs::GeV, 0.1 * cgs::PeV, 1000);
  auto units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;
  utils::OutputFile out("test_secondary_lis.txt");
  out << "# E [GeV] - H [] - He []\n";
  out << std::scientific;
  for (auto E : energyAxis) {
    out << E / cgs::GeV << " ";
    out << H_lis.get(E) / units << " ";
    out << He_lis.get(E) / units << " ";
    out << "\n";
  }
}

void printXsecs() {
  auto posXsecs = Orusa2022::PositronXsecs();
  auto xAxis = utils::LogAxis<double>(1e-4, 1, 1000);
  auto units = cgs::mbarn;
  auto channel = Orusa2022::Channel::Hep;
  utils::OutputFile out("test_secondary_xsecs.txt");
  out << "# x\n";
  out << std::scientific;
  for (auto x : xAxis) {
    out << x << " ";
    auto E = 10. * cgs::GeV;
    out << E * posXsecs.get(channel, E, x * E) / units << " ";
    E = 100. * cgs::GeV;
    out << E * posXsecs.get(channel, E, x * E) / units << " ";
    E = cgs::TeV;
    out << E * posXsecs.get(channel, E, x * E) / units << " ";
    E = 10. * cgs::TeV;
    out << E * posXsecs.get(channel, E, x * E) / units << " ";
    out << "\n";
  }
}

void printSourceTerm() {
  auto in = core::Input();
  in.print();
  auto sec = particle::SecondaryPositrons(in);
  auto energyAxis = utils::LogAxis<double>(1e-2 * cgs::GeV, 1e3 * cgs::GeV, 100);
  auto units = 1. / cgs::GeV / cgs::cm3 / cgs::sec;
  utils::OutputFile out("test_secondary_sourceterm.txt");
  out << "# E [GeV]\n";
  out << std::scientific;
  for (auto E : energyAxis) {
    out << E / cgs::GeV << " ";
    out << sec.Q_pp(E) / units << " ";
    out << sec.Q_Hep(E) / units << " ";
    out << "\n";
  }
}

void printIntensity() {
  auto in = core::Input();
  in.set_refEnergy(10. * cgs::GeV);
  in.set_D0_over_H(1.3e28 * cgs::cm2 / cgs::sec / cgs::kpc);  // At 10 GeV!
  in.set_delta(0.56);
  in.set_halosize(5. * cgs::kpc);
  in.set_galaxyRadius(20. * cgs::kpc);
  in.print();
  auto sec = particle::SecondaryPositrons(in);
  auto energyAxis = utils::LogAxis<double>(1e-1 * cgs::GeV, 1e4 * cgs::GeV, 100);
  auto units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;
  utils::OutputFile out("test_secondary_intensity_H5.txt");
  out << "# E [GeV]\n";
  out << std::scientific;
  for (auto E : energyAxis) {
    out << E / cgs::GeV << " ";
    out << sec.get(E) / units << " ";
    out << "\n";
  }
}

int main() {
  try {
    utils::startup_information();
    // printSourceTerm();
    printIntensity();
  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
  }
  return 0;
}