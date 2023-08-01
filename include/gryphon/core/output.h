#ifndef GRYPHON_CORE_OUTPUT_H
#define GRYPHON_CORE_OUTPUT_H

#include <memory>
#include <vector>

#include "gryphon/core/cgs.h"
#include "gryphon/core/input.h"
#include "gryphon/galaxy/galaxy.h"
#include "gryphon/particle/particle.h"

namespace gryphon {
namespace core {

class OutputManager {
 private:
  std::vector<double> m_E;
  std::vector<double> m_I;
  std::string m_filename;

 public:
  OutputManager(const Input& input);
  virtual ~OutputManager() = default;
  void dump() const;
  void compute(const std::shared_ptr<galaxy::Galaxy>& galaxy,
               const std::shared_ptr<particle::Particle>& particle);
};

}  // namespace core
}  // namespace gryphon

#endif