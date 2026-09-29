// Single-power-law injection spectrum for a few slopes, at the cutoff energy
// and efficiency of the configuration.
#include <array>
#include <cstdlib>
#include <iomanip>

#include "gryphon.h"

using namespace gryphon;

namespace {

constexpr std::array<double, 3> kSlopes{{1.5, 2.0, 2.3}};

}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "inspectInjectionSpectrum");
    if (!args) return EXIT_SUCCESS;

    auto input = utils::makeInput(*args);
    input.print();

    std::vector<injection::SinglePowerLawSpectrum> spectra;
    for (const auto slope : kSlopes) {
      input.set_injSlope(slope);
      spectra.emplace_back(input);
    }

    const core::RunOutput run(input, args->outdir);
    auto out = run.open("injectionspectrum",
                        "E [GeV] | E2Q(alpha=1.5) [erg] | E2Q(alpha=2.0) [erg] | "
                        "E2Q(alpha=2.3) [erg]");
    out << std::scientific << std::setprecision(6);

    for (const auto E : utils::LogAxis<double>(cgs::GeV, 1e6 * cgs::GeV, 1000)) {
      out << E / cgs::GeV;
      for (const auto& spectrum : spectra) {
        out << "\t" << (pow2(E) * spectrum.get(E)) / cgs::erg;
      }
      out << "\n";
    }

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
