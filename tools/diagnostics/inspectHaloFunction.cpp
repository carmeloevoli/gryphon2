// Halo function across the halo, for a range of diffusion lengths.
#include <array>
#include <cstdlib>
#include <iomanip>

#include "gryphon.h"

using namespace gryphon;

namespace {

// Diffusion lengths probed, in units of the halo half-height.
constexpr std::array<double, 5> kLambdaOverH{{0.1, 0.3, 0.5, 1.0, 10.0}};

}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "inspectHaloFunction");
    if (!args) return EXIT_SUCCESS;

    const auto input = utils::makeInput(*args);
    input.print();

    const core::RunOutput run(input, args->outdir);
    auto out = run.open("halofunction",
                        "z [kpc] | f(lambda=0.1H) | f(lambda=0.3H) | f(lambda=0.5H) | "
                        "f(lambda=1H) | f(lambda=10H)");
    out << std::scientific << std::setprecision(6);

    for (const auto z : utils::LinAxis<double>(-input.H(), input.H(), 1001)) {
      out << z / cgs::kpc;
      for (const auto ratio : kLambdaOverH) {
        out << "\t" << utils::halo_function(pow2(ratio * input.H()), input.H(), z, 0.);
      }
      out << "\n";
    }

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
