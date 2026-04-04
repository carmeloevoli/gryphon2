#include "gryphon/core/xsecs.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "gryphon/core/cgs.h"
#include "gryphon/utils/io.h"
#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace Orusa2022 {

namespace {

bool isSkippableLine(const std::string& line) {
  const auto first = line.find_first_not_of(" \t\r\n");
  return first == std::string::npos || line[first] == '#';
}

std::vector<double> parseColumns(const std::string& line) {
  std::istringstream stream(line);
  std::vector<double> values;
  double value = 0.;
  while (stream >> value) values.push_back(value);

  if (!stream.eof()) {
    throw std::runtime_error("failed to parse cross-section row: " + line);
  }

  return values;
}

}  // namespace

XSECS::XSECS(const std::string& filename) { read_filename(filename); }

void XSECS::read_filename(const std::string& filename) {
  double TprojMin = 0.1 * cgs::GeV;
  double TprojMax = 1e6 * cgs::GeV;
  size_t TprojSize = 140;
  m_lgTprojAxis = utils::LinAxis(std::log(TprojMin), std::log(TprojMax), TprojSize);

  double TsecMin = 0.01 * cgs::GeV;
  double TsecMax = 1e4 * cgs::GeV;
  size_t TsecSize = 90;
  m_lgTsecAxis = utils::LinAxis(std::log(TsecMin), std::log(TsecMax), TsecSize);

  double units = cgs::mbarn / cgs::GeV;

  m_sigma_pp.clear();
  m_sigma_pHe.clear();
  m_sigma_Hep.clear();
  m_sigma_HeHe.clear();

  if (!utils::fileExists(filename)) {
    throw std::runtime_error("cross-section table cannot be found: " + filename);
  }

  std::string line;
  std::ifstream file(filename.c_str());
  if (!file.is_open()) {
    throw std::runtime_error("could not open cross-section table: " + filename);
  }

  while (getline(file, line)) {
    if (isSkippableLine(line)) continue;

    const auto columns = parseColumns(line);
    if (columns.size() <= 9) {
      throw std::runtime_error("cross-section row has too few columns in '" + filename + "'");
    }

    m_sigma_pp.push_back(columns[2] * units);
    m_sigma_pHe.push_back(columns[3] * units);
    m_sigma_Hep.push_back(columns[8] * units);
    m_sigma_HeHe.push_back(columns[9] * units);
  }

  const size_t expectedSize = TprojSize * TsecSize;
  if (m_sigma_pp.size() != expectedSize || m_sigma_pHe.size() != expectedSize ||
      m_sigma_Hep.size() != expectedSize || m_sigma_HeHe.size() != expectedSize) {
    throw std::runtime_error("cross-section table '" + filename + "' has " +
                             std::to_string(m_sigma_pp.size()) + " rows, expected " +
                             std::to_string(expectedSize));
  }
}

double XSECS::get(const Channel& ch, const double& T_proj, const double& T_lepton) const {
  double value = 0;
  if (T_lepton > T_proj) return 0.;

  const auto lgTproj = std::log(T_proj);
  const auto lgTsec = std::log(T_lepton);

  if (utils::isInside(lgTproj, m_lgTprojAxis) && utils::isInside(lgTsec, m_lgTsecAxis)) {
    if (ch == Channel::pp) {
      value = utils::interpolate2d(lgTproj, lgTsec, m_lgTprojAxis, m_lgTsecAxis, m_sigma_pp);
    } else if (ch == Channel::pHe) {
      value = utils::interpolate2d(lgTproj, lgTsec, m_lgTprojAxis, m_lgTsecAxis, m_sigma_pHe);
    } else if (ch == Channel::Hep) {
      value = utils::interpolate2d(lgTproj, lgTsec, m_lgTprojAxis, m_lgTsecAxis, m_sigma_Hep);
    } else if (ch == Channel::HeHe) {
      value = utils::interpolate2d(lgTproj, lgTsec, m_lgTprojAxis, m_lgTsecAxis, m_sigma_HeHe);
    }
  }
  return value;
}

}  // namespace Orusa2022
}  // namespace gryphon
