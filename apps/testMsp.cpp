#include "gryphon.h"

using namespace gryphon;

void galacticLuminosity(double P_0, double B_S) {
  auto I = 2. / 5. * cgs::pulsar_mass * pow2(cgs::pulsar_radius);
  LOGD << "I : " << I / (cgs::gram / cgs::cm2) << " g/cm2";
  auto Omega_0 = 2. * M_PI / P_0;
  auto tau_0 = 3. * pow3(cgs::c_light) * I / pow2(B_S) / pow6(cgs::pulsar_radius) / pow2(Omega_0);
  auto L = 0.5 * I * pow2(Omega_0) / tau_0;
  LOGD << "tau_0 : " << tau_0 / cgs::Gyr << " Gyr";
  LOGD << "L : " << L / (cgs::erg / cgs::sec) << " erg/s";
  LOGD << " MSP L : " << 3e5 * L / (cgs::erg / cgs::sec) << " erg/s";
  auto V = 2. * pow2(M_PI) * B_S * pow3(cgs::pulsar_radius) / pow2(cgs::c_light) / pow2(P_0);
  LOGD << "V : " << cgs::elementary_charge * V / cgs::TeV << " TeV";
  auto N = 3e5;
  auto t_H = 10. * cgs::Gyr;
  LOGD << "Rate : " << N / t_H * cgs::Myr << " Myr^-1";
  LOGD << " PWN L : " << 0.1 * (2. / 100. / cgs::year) * 1e48 * cgs::erg << " erg/s";
}

void printEstar() {
  auto in = core ::Input();
  in.print();
  auto particle = std::make_shared<particle::MSP>(in);
  auto units = cgs::GeV;
  {
    auto energyAxis = utils::LogAxis<double>(cgs::GeV, cgs::PeV, 10000);
    utils::OutputFile out("test_msp_Estar.txt");
    out << "# E [GeV] - E_star []\n";
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / cgs::GeV << "\t";
      out << particle->Estar(E, 1e1 * cgs::kyr) / units << "\t";
      out << particle->Estar(E, 1e2 * cgs::kyr) / units << "\t";
      out << particle->Estar(E, 1e3 * cgs::kyr) / units << "\t";
      out << particle->Estar(E, 1e4 * cgs::kyr) / units << "\t";
      out << particle->Estar(E, 1e5 * cgs::kyr) / units << "\t";
      out << "\n";
    }
  }
}

void printLambda() {
  auto in = core ::Input();
  in.print();
  auto particle = std::make_shared<particle::MSP>(in);
  {
    auto energyAxis = utils::LogAxis<double>(cgs::GeV, cgs::PeV, 10000);
    utils::OutputFile out("test_msp_lambda.txt");
    out << "# E [GeV] - E_star []\n";
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / cgs::GeV << "\t";
      out << std::sqrt(particle->lambda2(E, 1e10 * E)) / cgs::kpc << "\t";
      out << particle->tau(E, 1e10 * E) / cgs::Myr << "\t";
      out << "\n";
    }
  }
}

void printSingleSpectrum() {
  auto in = core ::Input();
  in.print();
  auto particle = std::make_shared<particle::MSP>(in);
  auto units = 1. / (1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr);
  const auto d = utils::Vector3d(cgs::kpc, 0., 0.);
  {
    auto energyAxis = utils::LogAxis<double>(0.1 * cgs::GeV, 10. * cgs::TeV, 10000);
    utils::OutputFile out("test_msp_spectrum.txt");
    out << "# E [GeV] - flux []\n";
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / cgs::GeV << "\t";
      out << particle->get(E, 1e2 * cgs::kyr, d) / units << "\t";
      out << particle->get(E, 1e3 * cgs::kyr, d) / units << "\t";
      out << particle->get(E, 1e4 * cgs::kyr, d) / units << "\t";
      out << particle->get(E, 1e5 * cgs::kyr, d) / units << "\t";
      out << particle->get(E, 1e6 * cgs::kyr, d) / units << "\t";
      out << particle->get(E, 1e7 * cgs::kyr, d) / units << "\t";
      out << particle->Q(E) << "\t";
      out << "\n";
    }
  }
}

void run(unsigned long int seed, std::string simName) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.set_simname(simName);
  in.set_simEmin(cgs::GeV);
  in.set_simEmax(10. * cgs::TeV);
  in.set_simEsize(4 * 16);
  in.set_rate(30. / cgs::Myr);
  in.set_maxtime(10. * cgs::Gyr);
  in.print();

  RandomNumberGenerator rng = utils::RNG<double>(in.seed);

  auto galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
  galaxy->generate(rng, false);

  particle::Particles particles;
  particles.reserve(galaxy->size());

  auto events = galaxy->get_events();
  for (auto& event : events) {
    auto particle = std::make_shared<particle::MSP>(in, event, rng);
    particles.emplace_back(particle);
  }

  auto output = std::make_shared<core::OutputManager>(in);
  output->compute(particles);
  output->dump();
}

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    if (argc != 2) throw std::runtime_error("Usage: ./run params.ini");
    utils::Timer timer("timer for main");
    galacticLuminosity(5 * cgs::msec, 1e8 * cgs::gauss);
    // printEstar();
    // printLambda();
    // printSingleSpectrum();
    run(atoi(argv[1]), "msp");

  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
  }
  return 0;
}