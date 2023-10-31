#ifndef GRYPHON_PARTICLE_SECONDARY_H
#define GRYPHON_PARTICLE_SECONDARY_H

#include "gryphon/core/input.h"
#include "gryphon/core/lis.h"
#include "gryphon/core/xsecs.h"
#include "gryphon/particle/particle.h"

namespace gryphon {
namespace particle {

class SecondaryPositrons final : public Particle {
 public:
  SecondaryPositrons(const core::Input& in);
  virtual ~SecondaryPositrons() = default;

  double Q(double E) const override;
  double Q_pp(double E) const;
  double Q_Hep(double E) const;
  double get(double E) const override;

 protected:
  double b(double E) const;
  double tau(double E, double E_s) const;
  double Estar(double E, double dt) const;
  double lambda2(double E, double E_s) const;

 protected:
  const double m_E0 = 10. * cgs::GeV;
  const double m_D0 = 0.2 * cgs::kpc2 / cgs::Myr;
  const double m_delta = 0.56;
  const double m_H = 10. * cgs::kpc;
  const double m_b0 = 0.3 * cgs::GeV / cgs::Myr;
  const double m_Rd = 20. * cgs::kpc;

  Boschini2017::LIS H_lis = Boschini2017::H_LIS();
  Boschini2017::LIS He_lis = Boschini2017::He_LIS();
  Orusa2022::XSECS posXsecs = Orusa2022::PositronXsecs();
};

}  // namespace particle
}  // namespace gryphon

#endif