#ifndef GRYPHON_PARTICLE_VARYEFFICIENCY_H
#define GRYPHON_PARTICLE_VARYEFFICIENCY_H

#include "gryphon/particle/particle.h"
#include "gryphon/utils/random.h"

namespace gryphon {
namespace particle {

class VaryEfficiencyParticle final : public Particle {
 public:
  VaryEfficiencyParticle(const core::Input& in, RandomNumberGenerator& rng);
  virtual ~VaryEfficiencyParticle() = default;

  double Q(double E) const override;
  double get(double E, double dt, utils::Vector3d pos) const override;

 protected:
  double source_normalization() const;

 protected:
  RandomNumberGenerator& m_rng;
  const double m_E0 = 1. * cgs::GeV;
  double m_Q0 = 0;
  double m_varEfficienvy = 0;
};

}  // namespace particle
}  // namespace gryphon

#endif  // GRYPHON_PARTICLE_VARYEFFICIENCY_H