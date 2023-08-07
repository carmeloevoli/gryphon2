#ifndef GRYPHON_PARTICLE_H
#define GRYPHON_PARTICLE_H

#include <memory>
#include <vector>

#include "gryphon/core/event.h"
#include "gryphon/core/pid.h"
#include "gryphon/utils/vector3.h"

namespace gryphon {
namespace particle {

class Particle {
 public:
  Particle(const core::PID& pid, const std::shared_ptr<core::Event>& event);
  virtual ~Particle() = default;

  virtual double Q(double E) const = 0;
  virtual double get(double E) const = 0;

  const double& age = m_age;

 protected:
  core::PID m_pid;
  double m_age;
  utils::Vector3d m_pos;
};

using Particles = std::vector<std::shared_ptr<Particle> >;

}  // namespace particle
}  // namespace gryphon

#endif  // SIMPROP_PARTICLESTACK_H