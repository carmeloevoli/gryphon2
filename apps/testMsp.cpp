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

void printSingleSpectrum(const size_t id, const double distance, const double age,
                         const double luminosity) {
  auto in = core ::Input();
  // in.print();
  auto particle = std::make_shared<particle::MSP>(in);
  auto units = 1. / (1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr);
  const auto d = utils::Vector3d(distance, 0., 0.);
  {
    auto energyAxis = utils::LogAxis<double>(0.1 * cgs::GeV, 10. * cgs::TeV, 10000);
    utils::OutputFile out("test_msp_spectrum_" + std::to_string(id) + ".txt");
    out << "# E [GeV] - flux []\n";
    out << std::scientific;
    const auto factor = luminosity / particle->Luminosity;
    LOGD << "Scaling factor : " << factor;
    for (auto E : energyAxis) {
      out << E / cgs::GeV << "\t";
      out << factor * particle->get(E, age, d) / units << "\t";
      out << particle->Q(E) << "\t";
      out << "\n";
    }
  }
}

void printPropagation(const size_t id) {
  auto in = core ::Input();
  auto particle = std::make_shared<particle::MSP>(in);
  {
    auto energyAxis = utils::LogAxis<double>(cgs::GeV, cgs::PeV, 10000);
    utils::OutputFile out("test_msp_propagation_" + std::to_string(id) + ".txt");
    out << "# E [GeV] - tau [] - lambda []\n";
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / cgs::GeV << "\t";
      out << particle->tau(E, 1e4 * E) / cgs::Myr << "\t";
      out << std::sqrt(particle->lambda2(E, 1e4 * E)) / cgs::kpc << "\t";
      out << "\n";
    }
  }
}

// void run(unsigned long int seed, std::string simName) {
//   auto in = core ::Input();
//   in.set_seed(seed);
//   in.set_simname(simName);
//   in.set_simEmin(cgs::GeV);
//   in.set_simEmax(10. * cgs::TeV);
//   in.set_simEsize(64);
//   in.set_rate(30. / cgs::Myr);
//   in.set_maxtime(10. * cgs::Gyr);
//   in.print();

//   RandomNumberGenerator rng = utils::RNG<double>(in.seed);

//   auto galaxy = std::make_shared<galaxy::GalaxySteiman2010>(in);
//   galaxy->generate(rng, false);

//   particle::Particles particles;
//   particles.reserve(galaxy->size());

//   auto events = galaxy->get_events();
//   for (auto& event : events) {
//     auto particle = std::make_shared<particle::MSP>(in, event, rng);
//     particles.emplace_back(particle);
//   }

//   auto output = std::make_shared<core::OutputManager>(in);
//   output->compute(particles);
//   output->dump();
// }

struct Source {
  size_t id;
  double distance;  // from Earth
  double age;
  double luminosity;
};

std::vector<Source> readSources(const std::string& filename) {
  std::vector<Source> sources;
  std::ifstream file(filename);
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    // Skip comments and empty lines
    if (line.empty() || line[0] == '#') continue;
    Source source;
    if (iss >> source.id >> source.distance >> source.age >> source.luminosity) {
      sources.push_back(source);
    }
  }
  return sources;
}

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    if (argc != 2) throw std::runtime_error("Usage: ./run params.ini");
    utils::Timer timer("timer for main");
    std::string filename("top100_sources.txt");
    auto sources = readSources(filename);
    size_t count = 0;
    for (const auto& source : sources) {
      LOGD << "Processing source ID: " << source.id << ", Distance: " << source.distance
           << " kpc, Age: " << source.age << " Myr, Luminosity: " << source.luminosity << " erg/s";
      printSingleSpectrum(source.id, source.distance * cgs::kpc, source.age * cgs::Myr,
                          source.luminosity * cgs::erg / cgs::sec);
      if (count == 0) printPropagation(source.id);
      count++;
    }
    LOGD << "Processed " << count << " sources.";

    // galacticLuminosity(5 * cgs::msec, 1e8 * cgs::gauss);
    //   printEstar();
    //   printLambda();    //  run(atoi(argv[1]), "msp");

  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
  }
  return 0;
}