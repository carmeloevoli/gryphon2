#ifndef GRYPHON_DIFFUSION_H
#define GRYPHON_DIFFUSION_H

#include "gryphon/core/input.h"

namespace gryphon {
namespace core {

class DiffusionCoefficient {
 public:
  DiffusionCoefficient() {}
  DiffusionCoefficient(const Input& input);
  virtual ~DiffusionCoefficient() {}

  double get(const double& E) const;
  double escapeTimescale(const double& E) const;
  double lambda2(const double& E) const;

 protected:
  double m_H = 0;
  double m_D0 = 0;
  double m_E0 = 0;
  double m_delta = 0;
  double m_ddelta = 0;
  double m_s = 0;
  double m_Eb = 0;
};

}  // namespace core
}  // namespace gryphon

#endif  // GRYPHON_DIFFUSION_H
