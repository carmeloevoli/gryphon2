#include "gryphon/core/lis.h"

#include <cmath>

#include "gryphon/core/cgs.h"
#include "gryphon/utils/io.h"
#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace Boschini2017 {

LIS::LIS(core::PID pid, const std::string& filename) : Z(pid.get_Z()), A(pid.get_A()) {
  read_filename(filename);
}

double LIS::get(const double& T_n) const {
  double value = 0;
  double p = double(A) * std::sqrt(T_n * (T_n + 2. * cgs::proton_mass_c2));
  double logR_GV = std::log(p / double(Z) / cgs::GeV);
  if (logR_GV > logR_table.front() && logR_GV < logR_table.back()) {
    double logdJdR_GV = utils::interpolate(logR_GV, logR_table, logdJdR_table);
    double dJdR_GV = std::exp(logdJdR_GV);
    double dRdT = double(A) * (T_n + cgs::proton_mass_c2);
    dRdT /= double(Z) * std::sqrt(T_n * T_n + 2. * cgs::proton_mass_c2 * T_n);
    value = dRdT * dJdR_GV;
  }
  return value / (cgs::GeV * cgs::m2 * cgs::second);
}

void LIS::read_filename(const std::string& filename) {
  if (utils::fileExists(filename)) {
    std::fstream file_to_read(filename);
    const int max_num_of_char_in_a_line = 512;
    const int num_of_header_lines = 15;
    for (int i = 0; i < num_of_header_lines; ++i)
      file_to_read.ignore(max_num_of_char_in_a_line, '\n');
    double R_read, dJdR_read;
    while (!file_to_read.eof()) {
      file_to_read >> R_read >> dJdR_read;
      if (!file_to_read.eof()) {
        logR_table.push_back(std::log(R_read));
        logdJdR_table.push_back(std::log(dJdR_read));
      }
    }
  } else {
    throw std::runtime_error("file for reading LIS cannot be found.");
  }
}

}  // namespace Boschini2017
}  // namespace gryphon