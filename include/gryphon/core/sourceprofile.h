#ifndef GRYPHON_CORE_SOURCEPROFILE_H
#define GRYPHON_CORE_SOURCEPROFILE_H

#include "gryphon/core/input.h"

namespace gryphon {
namespace core {

class SourceProfile {
 public:
  SourceProfile(const Input& input);
  virtual ~SourceProfile() {}

  double get(const double& r) const;
  //   double escapeTimescale(const double& E) const;
  //   double lambda2(const double& E) const;

 protected:
  double m_a;
  double m_b;
  double m_R1;
  double m_Rsun;
};

}  // namespace core
}  // namespace gryphon

#endif  // GRYPHON_CORE_SOURCEPROFILE_H
