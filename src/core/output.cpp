#include "gryphon/core/output.h"

#include "gryphon/utils/io.h"
#include "gryphon/utils/logging.h"
#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace core {

std::string make_filename(std::string initfilename, int id, std::string ext) {
  std::string str = initfilename;
  str += "_";
  str += std::to_string(id);
  str += ".";
  str += ext;
  return str;
}

OutputManager::OutputManager(const Input& input) {
  m_E = utils::LogAxis<double>(input.E_min, input.E_max, input.E_size);
  m_filename = make_filename(input.simname, input.seed, "txt");
}

void OutputManager::dump() const {
  utils::OutputFile out(m_filename);
  out << "# E [GeV]\n";
  out << std::scientific;
  for (auto E : m_E) {
    out << E / cgs::GeV << "\t";
    // out << D.get(E) / (cgs::cm2 / cgs::sec) << "\t";
    out << "\n";
  }
}

//   std::ofstream outfile(fluxFilename.c_str());
//   if (outfile.is_open()) {
//     outfile << std::scientific << std::setprecision(4);
//     outfile << "# E [GeV] - flux [1 / GeV m^2 s sr]\n";
//     for (size_t i = 0; i < particle->get_size(); ++i) {
//       outfile << particle->get_E(i) / cgs::GeV << " ";
//       outfile << particle->get_total_flux(i) / units << " ";
//       outfile << particle->get_max_flux(i) / units << " ";
//       outfile << "\n";
//     }

}  // namespace core
}  // namespace gryphon