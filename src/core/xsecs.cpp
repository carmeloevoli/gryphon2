#include "gryphon/core/xsecs.h"

#include <cmath>

#include "gryphon/core/cgs.h"
#include "gryphon/utils/io.h"
#include "gryphon/utils/numeric.h"

namespace gryphon {
namespace Orusa2022 {

#define NHEADERLINES 42

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

  if (utils::fileExists(filename)) {
    std::string line;
    std::ifstream file(filename.c_str());
    while (getline(file, line)) {
      if (line.at(0) != '#') {
        auto s = utils::split(line, "             ");
        m_sigma_pp.push_back(std::stod(s[2]) * units);
        m_sigma_pHe.push_back(std::stod(s[3]) * units);
        m_sigma_Hep.push_back(std::stod(s[8]) * units);
        m_sigma_HeHe.push_back(std::stod(s[9]) * units);
      }
    }
  }
  assert(m_sigma_pp.size() == TprojSize * TsecSize);
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