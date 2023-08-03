#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
    auto in = core ::Input();
    in.print();
    auto zAxis = utils::LinAxis<double>(-in.H, in.H, 1001);
    utils::OutputFile out("test_halo.txt");
    out << "# z [kpc] - D [cm2/s]\n";
    out << std::scientific;
    for (auto z : zAxis) {
      out << z / cgs::kpc << "\t";
      out << utils::halo_function(pow2(0.1 * in.H), in.H, z, 0.) << "\t";
      out << utils::halo_function(pow2(0.5 * in.H), in.H, z, 0.) << "\t";
      out << utils::halo_function(pow2(in.H), in.H, z, 0.) << "\t";
      out << utils::halo_function(pow2(3. * in.H), in.H, z, 0.) << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return EXIT_SUCCESS;
}