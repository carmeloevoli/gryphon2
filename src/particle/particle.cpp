#include "gryphon/particle/particle.h"

namespace gryphon {
namespace particle {

Particle::Particle(const core::Input& in) : m_pid(in.pid), m_D(in) {
  m_alpha = in.injSlope;
  m_epsilon = in.injEfficiency;
  m_Emax = in.injEmax;
  m_H = in.H;
}

}  // namespace particle
}  // namespace gryphon