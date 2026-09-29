#include "gryphon/core/observables.h"

#include <iomanip>

#include "gryphon/core/cgs.h"

namespace gryphon {
namespace core {

void dumpFluxAndDipole(const RunOutput& run, const CosmicRays& cosmicRays) {
  const double flux_units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;
  auto out = run.open(
      "flux",
      "E [GeV] | I [GeV^-1 m^-2 s^-1 sr^-1] | dipole | dipole_x | dipole_y | dipole_z");
  out << std::scientific << std::setprecision(8);

  const auto& energy = cosmicRays.get_energyAxis();
  const auto& flux = cosmicRays.get_flux();
  const auto& dipole = cosmicRays.get_dipole();
  for (size_t i = 0; i < energy.size(); ++i) {
    out << energy[i] / cgs::GeV << "\t";
    out << flux[i] / flux_units << "\t";
    out << dipole[i].getModule() << "\t";
    out << dipole[i].x << "\t";
    out << dipole[i].y << "\t";
    out << dipole[i].z << "\n";
  }
}

}  // namespace core
}  // namespace gryphon
