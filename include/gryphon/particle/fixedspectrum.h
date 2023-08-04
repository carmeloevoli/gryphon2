#ifndef GRYPHON_PARTICLE_FIXEDSPECTRUM_H
#define GRYPHON_PARTICLE_FIXEDSPECTRUM_H

#include "gryphon/particle/particle.h"

namespace gryphon {
namespace particle {

class FixedSpectrumParticle final : public Particle {
 public:
  FixedSpectrumParticle(const core::Input& in);
  virtual ~FixedSpectrumParticle() = default;

  double Q(double E) const override;
  double get(double E, double dt, utils::Vector3d pos) const override;

 protected:
  double source_normalization() const;

 protected:
  const double m_E0 = 1. * cgs::GeV;
  double m_Q0 = 0;
};

}  // namespace particle
}  // namespace gryphon

#endif  // GRYPHON_PARTICLE_FIXEDSPECTRUM_H