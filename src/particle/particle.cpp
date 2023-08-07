#include "gryphon/particle/particle.h"

namespace gryphon {
namespace particle {

Particle::Particle(const core::PID& pid, const std::shared_ptr<core::Event>& event)
    : m_pid(pid), m_age(event->age), m_pos(event->pos) {}

}  // namespace particle
}  // namespace gryphon