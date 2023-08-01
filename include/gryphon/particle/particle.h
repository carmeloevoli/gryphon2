#ifndef GRYPHON_PARTICLE_H
#define GRYPHON_PARTICLE_H

#include "gryphon/core/diffusion.h"
#include "gryphon/core/input.h"
#include "gryphon/core/pid.h"
#include "gryphon/utils/vector3.h"

namespace gryphon {
namespace particle {

class Particle {
 protected:
  PID::PID m_pid;
  core::DiffusionCoefficient m_D;

 public:
  Particle(const core::Input& in) : m_pid(in.pid), m_D(in) {}
  virtual ~Particle() = default;

  virtual double get(double E, double dt, utils::Vector3d pos) const = 0;
};

}  // namespace particle
}  // namespace gryphon

#endif  // SIMPROP_PARTICLESTACK_H