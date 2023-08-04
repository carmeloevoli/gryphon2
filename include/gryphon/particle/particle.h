#ifndef GRYPHON_PARTICLE_H
#define GRYPHON_PARTICLE_H

#include "gryphon/core/diffusion.h"
#include "gryphon/core/input.h"
#include "gryphon/core/pid.h"
#include "gryphon/utils/vector3.h"

namespace gryphon {
namespace particle {

class Particle {
 public:
  Particle(const core::Input& in);
  virtual ~Particle() = default;

  virtual double Q(double E) const = 0;
  virtual double get(double E, double dt, utils::Vector3d pos) const = 0;

 protected:
  PID::PID m_pid;
  core::DiffusionCoefficient m_D;
  double m_alpha;
  double m_epsilon;
  double m_Emax;
  double m_H;
};

}  // namespace particle
}  // namespace gryphon

#endif  // SIMPROP_PARTICLESTACK_H