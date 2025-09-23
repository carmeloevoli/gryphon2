#include "gryphon.h"

using namespace gryphon;

void calc(double P_0, double B_S) {
  auto I = 2. / 5. * cgs::pulsar_mass * pow2(cgs::pulsar_radius);
  LOGD << "I : " << I / (cgs::gram / cgs::cm2) << " g/cm2";
  auto Omega_0 = 2. * M_PI / P_0;
  LOGD << "Omega_0 : " << Omega_0 / (1. / cgs::second) << " 1/s";
  auto V = B_S * pow3(cgs::pulsar_radius) * pow2(Omega_0) / pow2(cgs::c_light);
  LOGD << "V : " << cgs::elementary_charge * V / cgs::eV << " eV";
  auto N_dot =
      B_S * pow3(cgs::pulsar_radius) * pow2(Omega_0) / (4. * cgs::elementary_charge * cgs::c_light);
  LOGD << "N_dot : " << N_dot / (1. / cgs::second) << " 1/s";
  auto E = 1e20 * cgs::eV;
  auto N_protons = 3. * I * pow2(cgs::c_light) / 4. / cgs::elementary_charge / B_S /
                   pow3(cgs::pulsar_radius) / E;
  LOGD << "N_protons : " << N_protons * cgs::GeV << " protons per pulsar (GeV-1)";
}

void compute_flux() {
  // Pulsar parameters
  const auto I = 2. / 5. * cgs::pulsar_mass * pow2(cgs::pulsar_radius);
  const auto pulsar_rate = 1. / (50. * cgs::year);  // One pulsar per 50 years
  const auto B_S = std::pow(10., 12.) * cgs::gauss;
  const auto R_d = 10. * cgs::kpc;    // Radius of the Galactic disc
  const auto E_max = 2e15 * cgs::eV;  // Maximum energy of emitted particles
  const auto P_0 = 80. * cgs::msec;   // Initial spin period of the pulsar

  // Diagnostics
  calc(P_0, B_S);

  auto factor = 3. * I * pow2(cgs::c_light) * pulsar_rate;
  factor /= 8. * M_PI * pow2(R_d) * cgs::elementary_charge * B_S * pow3(cgs::pulsar_radius) * E_max;

  // Compute the spectrum of emitted particles
  auto E = utils::LogAxis<double>(cgs::GeV, 10 * cgs::PeV, 1000);
  auto flux = std::vector<double>(E.size(), 0.0);
  auto units = 1. / (cgs::GeV * cgs::m2 * cgs::second * cgs::sr);
  for (size_t i = 0; i < E.size(); ++i) {
    flux[i] = factor;
    const auto D_H = 0.42 * std::pow(E[i] / cgs::TeV, 0.36) * cgs::kpc / cgs::Myr;
    flux[i] /= D_H;
    flux[i] *= (E_max / E[i]) * std::exp(-E[i] / E_max);
    flux[i] *= cgs::c_light / (4. * M_PI);
  }
  std::ofstream outfile("output/hpwne_model.txt");
  outfile << "# E [GeV] / flux [GeV-1 m-2 s-1 sr-1]" << std::endl;
  outfile << std::scientific << std::setprecision(3);
  for (size_t i = 0; i < E.size(); ++i) {
    outfile << E[i] / cgs::GeV << " " << flux[i] / units << std::endl;
  }
  outfile.close();
}

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    if (argc != 2) throw std::runtime_error("Usage: ./run params.ini");
    utils::Timer timer("timer for main");

    // Run simulation
    compute_flux();

    // auto in = core::Input(argv[1]);
    // run(69, "HP");

  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
  }
  return 0;
}