#ifndef GRYPHON_PARTICLESTACK_SNRBUILDER_H
#define GRYPHON_PARTICLESTACK_SNRBUILDER_H

#include "gryphon/particle/particle.h"

namespace gryphon {
namespace particle {

class FixedSpectrumParticle final : public Particle {
 protected:
  core::DiffusionCoefficient m_D;
  const double m_E0 = 1. * cgs::GeV;
  double m_alpha;
  double m_epsilon;
  double m_Emax;
  double m_Q0;

 public:
  FixedSpectrumParticle(const core::Input& in) : Particle(in) {
    m_alpha = in.injSlope;
    m_epsilon = in.injEfficiency;
    m_Emax = in.injEmax;
    m_Q0 = compute_normalization();
    m_D = core::DiffusionCoefficient(in);
  }

  virtual ~FixedSpectrumParticle() = default;

  double Q(double E) const override;
  double get(double E, double dt, utils::Vector3d pos) const override;

 protected:
  double compute_normalization() const;
};

}  // namespace particle

}  // namespace gryphon

#endif  // GRYPHON_PARTICLESTACK_SNRBUILDER_H