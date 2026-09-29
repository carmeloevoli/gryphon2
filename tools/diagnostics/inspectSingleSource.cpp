// Everything about a single source at 1 kpc: its injection spectrum, the flux
// it produces at a range of ages, the time evolution at fixed energies, and the
// transport timescales.
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <vector>

#include "gryphon.h"

using namespace gryphon;

namespace {

const utils::Vector3d kSourcePosition(cgs::kpc, 0., 0.);
constexpr double kFluxUnits = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;

void dumpInjectionSpectrum(const injection::InjectionSpectrum& spectrum,
                           const core::RunOutput& run) {
  auto out = run.open("injectionspectrum", "E [GeV] | Q(E) [GeV^-1]");
  out << std::scientific << std::setprecision(6);

  const auto units = 1. / cgs::GeV;
  for (const auto E : utils::LogAxis<double>(cgs::GeV, cgs::TeV, 100)) {
    out << E / cgs::GeV << "\t" << spectrum.get(E) / units << "\n";
  }
}

void dumpFluxAtAges(const core::Input& in, const std::shared_ptr<const kernel::GreenKernel>& greenKernel,
                    RandomNumberGenerator& rng, const core::RunOutput& run) {
  const std::vector<double> ages = utils::LogAxis<double>(0.01 * cgs::Myr, 1. * cgs::Myr, 10);

  for (size_t i = 0; i < ages.size(); ++i) {
    core::Events events;
    events.emplace_back(std::make_shared<core::Event>(ages[i], kSourcePosition));

    auto eventSpectra = injection::makeInjectionSpectra(in, events, rng);
    core::CosmicRays cr(in, greenKernel, std::move(eventSpectra), events);
    cr.run();

    std::ostringstream columns;
    columns << "E [GeV] | I(age=" << ages[i] / cgs::Myr << " Myr) [GeV^-1 m^-2 s^-1 sr^-1]";

    auto out = run.open("age" + std::to_string(i), columns.str());
    out << std::scientific << std::setprecision(6);

    const auto& E = cr.get_energyAxis();
    const auto& I = cr.get_flux();
    for (size_t j = 0; j < E.size(); ++j) {
      out << E[j] / cgs::GeV << "\t" << I[j] / kFluxUnits << "\n";
    }
  }
}

void dumpTimeScan(const std::shared_ptr<const kernel::GreenKernel>& greenKernel,
                  const std::shared_ptr<const injection::InjectionSpectrum>& spectrum,
                  const core::RunOutput& run) {
  const std::vector<double> energies = {1e1 * cgs::GeV, 1e2 * cgs::GeV, cgs::TeV};

  std::ostringstream columns;
  columns << "t [Myr]";
  for (const auto E : energies) {
    columns << " | I(E=" << E / cgs::GeV << " GeV) [GeV^-1 m^-2 s^-1 sr^-1]";
  }

  auto out = run.open("timescan", columns.str());
  out << std::scientific << std::setprecision(6);

  for (const auto t : utils::LogAxis<double>(1e-3 * cgs::Myr, 10. * cgs::Myr, 1000)) {
    out << t / cgs::Myr;
    for (const auto E : energies) {
      const double flux = greenKernel->flux(
          E, t, kSourcePosition, [spectrum](double Eprime) { return spectrum->get(Eprime); });
      out << "\t" << flux / kFluxUnits;
    }
    out << "\n";
  }
}

void dumpTimescales(const core::Input& in,
                    const std::shared_ptr<const kernel::GreenKernel>& greenKernel,
                    const core::RunOutput& run) {
  auto out = run.open("timescales", "E [GeV] | t_diff [Myr] | t_loss [Myr]");
  out << std::scientific << std::setprecision(6);

  for (const auto E : utils::LogAxis<double>(in.simEmin(), in.simEmax(), in.simEsize())) {
    out << E / cgs::GeV << "\t";
    out << greenKernel->diffusionTimescale(E) / cgs::Myr << "\t";
    out << greenKernel->energyLossTimescale(E) / cgs::Myr << "\n";
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "inspectSingleSource");
    if (!args) return EXIT_SUCCESS;

    const auto input = utils::makeInput(*args);
    input.print();

    RandomNumberGenerator rng(input.seed());
    auto greenKernel = kernel::makeGreenKernel(input);
    auto injectionSpectrum = injection::makeInjectionSpectrum(input, rng);

    const core::RunOutput run(input, args->outdir);
    dumpInjectionSpectrum(*injectionSpectrum, run);
    dumpFluxAtAges(input, greenKernel, rng, run);
    dumpTimeScan(greenKernel, injectionSpectrum, run);
    dumpTimescales(input, greenKernel, run);

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
