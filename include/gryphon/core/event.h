#ifndef GRYPHON_CORE_EVENT_H
#define GRYPHON_CORE_EVENT_H

#include <memory>
#include <vector>

#include "gryphon/utils/vector3.h"

namespace gryphon {
namespace core {

class Event {
 public:
  Event() {}
  Event(double age, utils::Vector3d posToSun) : m_age(age), m_pos(posToSun) {}
  virtual ~Event() = default;

 public:
  const double& age = m_age;
  const utils::Vector3d& pos = m_pos;

 private:
  double m_age;
  utils::Vector3d m_pos;
};

using Events = std::vector<std::shared_ptr<Event> >;

}  // namespace core
}  // namespace gryphon

#endif  // GRYPHON_CORE_EVENT_H