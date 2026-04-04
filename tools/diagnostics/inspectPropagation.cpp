#include <cassert>

#include "gryphon.h"

using namespace gryphon;

namespace {

void PureDiffusionDiagnostics(const core::Input& in) {
  auto kernel = std::make_shared<kernel::PureDiffusionKernel>(in);
  auto energyAxis = utils::LogAxis<double>(1e3 * cgs::GeV, 1e6 * cgs::GeV, 100);
  utils::OutputFile out("inspect_pure_diffusion.txt");
  out << "# E [GeV] - D [cm2/s] - t_diff [Myr] - n_sources\n";
  out << std::scientific;
  for (auto E : energyAxis) {
    const auto D = kernel->D(E);
    const auto t_diff = pow2(in.H()) / 2. / D;
    const auto n_sources = pow2(in.H() / in.R_g()) * in.sn_rate() * t_diff;
    out << E / cgs::GeV << "\t";
    out << D / (cgs::cm2 / cgs::sec) << "\t";
    out << t_diff / cgs::Myr << "\t";
    out << n_sources << "\t";
    out << "\n";
  }
}

void DiffusionLossesDiagnostics(const core::Input& in) {
  auto kernel = std::make_shared<kernel::DiffusionLossesKernel>(in);
  auto energyAxis = utils::LogAxis<double>(1e3 * cgs::GeV, 1e6 * cgs::GeV, 100);
  utils::OutputFile out("inspect_diffusion_losses.txt");
  out << "# E [GeV] - D [cm2/s]\n";
  out << std::scientific;
  for (auto E : energyAxis) {
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
    out << n_sources << "\t";
    out << "\n";
  }
}

}  // namespace

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    // in.set_refEnergy(cgs::GeV);
    // in.set_D0_over_H(0.35e28 * cgs::cm2 / cgs::sec / cgs::kpc);
    // in.set_delta(0.56);
    // in.set_halosize(5. * cgs::kpc);
    // in.set_galaxyRadius(20. * cgs::kpc);
    // in.set_Bfield(5. * cgs::microgauss);
    in.print();

    PureDiffusionDiagnostics(in);
    DiffusionLossesDiagnostics(in);
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}