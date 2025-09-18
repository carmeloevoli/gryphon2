#ifndef GRYPHON_CORE_PROPAGATION_H
#define GRYPHON_CORE_PROPAGATION_H

#include "gryphon/core/cgs.h"
#include "gryphon/core/input.h"

namespace gryphon {
namespace core {

class Propagation {
 public:
  Propagation() {}
  Propagation(const Input& in);
  virtual ~Propagation() = default;

  inline double D(double E) const { return m_D0 * std::pow(E / m_E0, m_delta); }
  inline double b(double E) const { return m_b0 * pow2(E / m_E0); }

 public:
  double tau(double E, double E_s) const;
  double Estar(double E, double dt) const;
  double lambda2(double E, double E_s) const;
  const double& H = m_H;

 protected:
  double m_E0 = cgs::GeV;
  double m_D0 = 1e28 * cgs::cm2 / cgs::sec;
  double m_H = cgs::kpc;
  double m_delta = 0.56;
  double m_b0 = 0.2 * cgs::GeV / cgs::Myr;
};

}  // namespace core
}  // namespace gryphon

#endif