#ifndef GRYPHON_PARTICLESTACK_SNRBUILDER_H
#define GRYPHON_PARTICLESTACK_SNRBUILDER_H

#include "gryphon/particle/particle.h"

namespace gryphon {
namespace particle {

class FixedSpectrumParticle final : public Particle {
 public:
  FixedSpectrumParticle(const core::Input& in) : Particle(in) {}
  virtual ~FixedSpectrumParticle() = default;

  double get(double E, double dt, utils::Vector3d pos) const override;
};

}  // namespace particle

}  // namespace gryphon

#endif  // GRYPHON_PARTICLESTACK_SNRBUILDER_H