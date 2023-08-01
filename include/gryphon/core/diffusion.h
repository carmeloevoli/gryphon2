#ifndef GRYPHON_DIFFUSION_H
#define GRYPHON_DIFFUSION_H

#include "gryphon/core/input.h"

namespace gryphon {
namespace core {

class DiffusionCoefficient {
 public:
  DiffusionCoefficient(const Input& input);
  virtual ~DiffusionCoefficient() {}

  double get(const double& E) const;
  double escapeTimescale(const double& E) const;
  double lambda2(const double& E) const;

 protected:
  double m_H;
  double m_D0;
  double m_E0;
  double m_delta;
  double m_ddelta;
  double m_s;
  double m_Eb;
};

}  // namespace core
}  // namespace gryphon

#endif  // GRYPHON_DIFFUSION_H
