#include "gryphon.h"

using namespace gryphon;

void calc(double P_0, double B_S) {
  auto I = 2. / 5. * cgs::pulsar_mass * pow2(cgs::pulsar_radius);
  LOGD << "I : " << I / (cgs::gram / cgs::cm2) << " g/cm2";
  auto Omega_0 = 2. * M_PI / P_0;
  auto tau_0 = 3. * pow3(cgs::c_light) * I / pow2(B_S) / pow6(cgs::pulsar_radius) / pow2(Omega_0);
  LOGD << "tau_0 : " << tau_0 / cgs::kyr << " kyr";
  auto L = 0.5 * I * pow2(Omega_0) / tau_0;
  LOGD << "L : " << L / (cgs::erg / cgs::sec) << " erg/s";
  auto V = 2. * pow2(M_PI) * B_S * pow3(cgs::pulsar_radius) / pow2(cgs::c_light) / pow2(P_0);
  LOGD << "V : " << cgs::elementary_charge * V / cgs::TeV << " TeV";
}

void run(unsigned long int seed, std::string simName) {
  auto in = core ::Input();
  in.set_seed(seed);
  in.set_simname(simName);
  in.set_simEmin(10. * cgs::GeV);
  in.set_simEmax(10. * cgs::PeV);
  in.set_simEsize(6 * 16);
  in.set_rate(0.03 / cgs::year);
  in.set_maxtime(cgs::Myr);
  in.print();

  RandomNumberGenerator rng = utils::RNG<double>(in.seed);

  auto galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
  galaxy->generate(rng, false);

  particle::Particles particles;
  particles.reserve(galaxy->size());

  //   auto events = galaxy->get_events();
  //   for (auto& event : events) {
  //     auto particle = std::make_shared<particle::MSP>(in, event, rng);
  //     particles.emplace_back(particle);
  //   }

  //   auto output = std::make_shared<core::OutputManager>(in);
  //   output->compute(particles);
  //   output->dump();
}

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    if (argc != 2) throw std::runtime_error("Usage: ./run params.ini");
    utils::Timer timer("timer for main");

    // Diagnostic for hadronic pulsars
    double P_0 = 0.1 * cgs::second;
    double B_S = std::pow(10., 12.65) * cgs::gauss;
    calc(P_0, B_S);

    // Run simulation
    auto in = core::Input(argv[1]);
    run(69, "HP");

  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
  }
  return 0;
}