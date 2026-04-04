#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.set_injectionModel(InjectionModel::SinglePowerLaw);
    in.set_transportModel(TransportModel::PureDiffusion);
    in.print();

    RandomNumberGenerator rng = utils::RNG<double>(in.seed());
    auto kernel = kernel::makeGreenKernel(in);
    auto injectionSpectrum = injection::makeInjectionSpectrum(in, rng);

    const auto d = utils::Vector3d(cgs::kpc, 0., 0.);
    const double units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;

    // 1) Spectrum for single-event solutions at different event ages.
    const std::vector<double> ages = utils::LogAxis<double>(0.01 * cgs::Myr, 0.1 * cgs::Myr, 10);

    for (size_t i = 0; i < ages.size(); ++i) {
      core::Events events;
      events.emplace_back(std::make_shared<core::Event>(ages[i], d));

      core::CosmicRays cr(in, kernel, injectionSpectrum, events);
      cr.run();

      utils::OutputFile out("test_single_age_" + std::to_string(i) + ".txt");
      out << "# age [Myr] = " << ages[i] / cgs::Myr << "\n";
      out << "# E [GeV] - I [GeV-1 m-2 sec-1 sr-1]\n";
      out << std::scientific;
      const auto& E = cr.get_energyAxis();
      const auto& I = cr.get_flux();
      for (size_t j = 0; j < E.size(); ++j) {
        out << E[j] / cgs::GeV << "\t";
        out << I[j] / units << "\n";
      }
    }

    // 2) Time evolution at fixed energies using the kernel directly.
    const std::vector<double> energies = {0.1 * cgs::TeV, cgs::TeV, 10. * cgs::TeV, 1e2 * cgs::TeV};
    const auto timeAxis = utils::LogAxis<double>(1e-3 * cgs::Myr, 10. * cgs::Myr, 1000);

    utils::OutputFile outTime("test_single_time_scan.txt");
    outTime << "# t [Myr]";
    for (auto E : energies) outTime << "\tI(E=" << E / cgs::GeV << " GeV)";
    outTime << "\n";
    outTime << std::scientific;

    for (auto t : timeAxis) {
      outTime << t / cgs::Myr;
      for (auto E : energies) {
        const double flux = kernel->flux(
            E, t, d, [injectionSpectrum](double Eprime) { return injectionSpectrum->get(Eprime); });
        outTime << "\t" << flux / units;
      }
      outTime << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}