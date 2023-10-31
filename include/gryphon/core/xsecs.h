#ifndef GRYPHON_CORE_XSECS_H
#define GRYPHON_CORE_XSECS_H

#include <string>
#include <vector>

namespace gryphon {
namespace Orusa2022 {

enum class Channel { pp, pHe, Hep, HeHe };

class XSECS {
 public:
  XSECS(const std::string &filename);
  virtual ~XSECS() = default;

  double get(const Channel &ch, const double &T_proj, const double &T_lepton) const;

 protected:
  void read_filename(const std::string &filename);

 private:
  std::vector<double> m_lgTprojAxis;
  std::vector<double> m_lgTsecAxis;
  std::vector<double> m_sigma_pp;
  std::vector<double> m_sigma_Hep;
  std::vector<double> m_sigma_pHe;
  std::vector<double> m_sigma_HeHe;
};

class PositronXsecs final : public XSECS {
 public:
  PositronXsecs() : XSECS("data/supplementary_table_positrons_best_fit.txt") {}
};

class ElectronXsecs final : public XSECS {
 public:
  ElectronXsecs() : XSECS("data/supplementary_table_electrons_best_fit.txt") {}
};

}  // namespace Orusa2022
}  // namespace gryphon

#endif