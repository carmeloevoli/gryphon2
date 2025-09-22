#ifndef GRYPHON_PARTICLE_HPWNE_H
#define GRYPHON_PARTICLE_HPWNE_H

#include "gryphon/core/event.h"
#include "gryphon/core/input.h"
#include "gryphon/core/propagation.h"
#include "gryphon/particle/particle.h"
#include "gryphon/utils/random.h"

namespace gryphon {
namespace particle {

class HPWNE final : public Particle {
 public:
  HPWNE(const core::Input& in, const std::shared_ptr<core::Event>& event,
        RandomNumberGenerator& rng);
  HPWNE(const core::Input& in);
  virtual ~HPWNE() = default;

  double Q(double E) const override;
  double get(double E, double dt, utils::Vector3d pos) const;
  double get(double E) const override;

 protected:
  const double m_E0 = 1. * cgs::GeV;
  double m_H = 0;

 protected:
  core::Propagation m_prop;
};

}  // namespace particle
}  // namespace gryphon

#endif  // GRYPHON_PARTICLE_HPWNE_H