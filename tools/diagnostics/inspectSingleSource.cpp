#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.set_simEmin(10. * cgs::GeV);
    in.set_simEmax(100. * cgs::TeV);
    in.set_simEsize(1000);
    in.set_pwnAlpha1(1.8);
    in.set_pwnAlpha2(1.8);
    in.set_pwnP0(30. * cgs::msec);
    in.set_efficiency(0.27);
    in.set_pwnEbreak(100. * cgs::GeV);
    in.set_injectionModel(InjectionModel::PWN);
    in.set_transportModel(TransportModel::DiffusionLosses);
    in.print();

    RandomNumberGenerator rng(in.seed());
    auto kernel = kernel::makeGreenKernel(in);
    auto injectionSpectrum = injection::makeInjectionSpectrum(in, rng);

    // 1) Injection spectrum for a single source.
    {
      utils::OutputFile out("inspect_single_source_injection_spectrum.txt");
      out << "# E [GeV] - Q [GeV-1]\n";
      out << std::scientific;
      const auto units = 1. / cgs::GeV;
      const auto energyAxis = utils::LogAxis<double>(cgs::GeV, cgs::TeV, 100);
      for (auto E : energyAxis) {
        out << E / cgs::GeV << "\t" << injectionSpectrum->get(E) / units << "\n";
      }
    }

    // 2) Spectrum for single-event solutions at different event ages.
    {
      const std::vector<double> ages = utils::LogAxis<double>(0.01 * cgs::Myr, 1 * cgs::Myr, 10);
      const auto d = utils::Vector3d(cgs::kpc, 0., 0.);
      const double units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;

      for (size_t i = 0; i < ages.size(); ++i) {
        core::Events events;
        events.emplace_back(std::make_shared<core::Event>(ages[i], d));

        auto eventSpectra = injection::makeInjectionSpectra(in, events, rng);
        core::CosmicRays cr(in, kernel, std::move(eventSpectra), events);
        cr.run();

        utils::OutputFile out("inspect_single_source_age_" + std::to_string(i) + ".txt");
        out << "# age [Myr] = " << ages[i] / cgs::Myr << "\n";
        out << "# E [GeV] - I [GeV-1 m-2 sec-1 sr-1]\n";
        out << std::scientific;
        const auto& E = cr.get_energyAxis();
        const auto& I = cr.get_flux();
        for (size_t j = 0; j < E.size(); ++j) {
          out << E[j] / cgs::GeV << "\t";
          out << I[j] / units << "\t";
          out << "\n";
        }
      }

      // 3) Time evolution at fixed energies using the kernel directly.
      {
        const std::vector<double> energies = {1e1 * cgs::GeV, 1e2 * cgs::GeV, cgs::TeV};
        const auto timeAxis = utils::LogAxis<double>(1e-3 * cgs::Myr, 10. * cgs::Myr, 1000);

        utils::OutputFile out("inspect_single_source_time_scan.txt");
        out << "# t [Myr]";
        for (auto E : energies) out << "\tI(E=" << E / cgs::GeV << " GeV)";
        out << "\n";
        out << std::scientific;

        for (auto t : timeAxis) {
          out << t / cgs::Myr;
          for (auto E : energies) {
            const double flux = kernel->flux(E, t, d, [injectionSpectrum](double Eprime) {
              return injectionSpectrum->get(Eprime);
            });
            out << "\t" << flux / units;
          }
          out << "\n";
        }
      }

      // 4) Timescales as a function of energy.
      {
        const auto energyAxis = utils::LogAxis<double>(in.simEmin(), in.simEmax(), in.simEsize());

        utils::OutputFile out("inspect_single_source_timescales.txt");
        out << "# E [GeV] - t_diff [Myr] - t_loss [Myr]\n";
        out << std::scientific;
        for (auto E : energyAxis) {
          const double t_diff = kernel->diffusionTimescale(E);
          const double t_loss = kernel->energyLossTimescale(E);
          out << E / cgs::GeV << "\t" << t_diff / cgs::Myr << "\t" << t_loss / cgs::Myr << "\n";
        }
      }
    }

    return EXIT_SUCCESS;
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
    return EXIT_FAILURE;
  }
}
