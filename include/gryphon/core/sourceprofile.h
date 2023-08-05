#ifndef GRYPHON_CORE_SOURCEPROFILE_H
#define GRYPHON_CORE_SOURCEPROFILE_H

#include "gryphon/core/input.h"
#include "gryphon/utils/lookupContainers.h"

namespace gryphon {
namespace core {

class SourceProfile {
 public:
  SourceProfile(const Input& input);
  virtual ~SourceProfile() = default;

  double get(const double& r) const;
  double pick_radius(const double& rnd) const;

 protected:
  double m_a;
  double m_b;
  double m_R1;
  double m_Rsun;
  double m_Rg;
  utils::LookupArray<20000> m_profileIntegral;

 protected:
  double integrateProfile(const double& r) const;
};

}  // namespace core
}  // namespace gryphon

#endif  // GRYPHON_CORE_SOURCEPROFILE_H
