// Transport diagnostics of a configuration: diffusion coefficient, escape and
// loss timescales, and the number of sources contributing at each energy.
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <vector>

#include "gryphon.h"

using namespace gryphon;

namespace {

std::vector<double> diagnosticEnergyAxis() {
  return utils::LogAxis<double>(1e3 * cgs::GeV, 1e6 * cgs::GeV, 100);
}

void dumpPureDiffusion(const core::Input& in, const core::RunOutput& run) {
  const auto kernel = std::make_shared<kernel::PureDiffusionKernel>(in);

  auto out = run.open("propagation",
                      "E [GeV] | D [cm2/s] | t_diff [Myr] | n_sources");
  out << std::scientific << std::setprecision(6);

  for (const auto E : diagnosticEnergyAxis()) {
    const auto D = kernel->D(E);
    const auto t_diff = pow2(in.H()) / 2. / D;
    const auto n_sources = pow2(in.H() / in.R_g()) * in.sn_rate() * t_diff;
    out << E / cgs::GeV << "\t";
    out << D / (cgs::cm2 / cgs::sec) << "\t";
    out << t_diff / cgs::Myr << "\t";
    out << n_sources << "\n";
  }
}

void dumpDiffusionLosses(const core::Input& in, const core::RunOutput& run) {
  const auto kernel = std::make_shared<kernel::DiffusionLossesKernel>(in);

  auto out = run.open("propagation",
                      "E [GeV] | D [cm2/s] | b [GeV/s] | t_diff [Myr] | t_loss [Myr] | "
                      "lambda [kpc] | n_sources");
  out << std::scientific << std::setprecision(6);

  for (const auto E : diagnosticEnergyAxis()) {
    const auto D = kernel->D(E);
    const auto b = kernel->b(E);
    const auto lambda2 = kernel->lambda2(E, 1e4 * E);
    const auto t_diff = pow2(in.H()) / 2. / D;
    const auto t_loss = E / b;
    const auto n_sources = pow2(in.H() / in.R_g()) * in.sn_rate() * std::min(t_diff, t_loss);
    out << E / cgs::GeV << "\t";
    out << D / (cgs::cm2 / cgs::sec) << "\t";
    out << b / (cgs::GeV / cgs::sec) << "\t";
    out << t_diff / cgs::Myr << "\t";
    out << t_loss / cgs::Myr << "\t";
    out << std::sqrt(lambda2) / cgs::kpc << "\t";
    out << n_sources << "\n";
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    utils::startup_information();
    const auto args = utils::parseCommandLine(argc, argv, "inspectPropagation");
    if (!args) return EXIT_SUCCESS;

    const auto input = utils::makeInput(*args);
    input.print();

    const core::RunOutput run(input, args->outdir);
    if (input.transportModel() == TransportModel::DiffusionLosses) {
      dumpDiffusionLosses(input, run);
    } else {
      dumpPureDiffusion(input, run);
    }

  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
