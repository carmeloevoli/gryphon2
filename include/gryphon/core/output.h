#ifndef GRYPHON_CORE_OUTPUT_H
#define GRYPHON_CORE_OUTPUT_H

#include <vector>

#include "gryphon/core/cgs.h"
#include "gryphon/core/input.h"

namespace gryphon {
namespace core {

class OutputManager {
 private:
  std::vector<double> m_E;
  std::string m_filename;
  const double units = 1. / cgs::GeV / cgs::m2 / cgs::sec / cgs::sr;

 public:
  OutputManager(const Input& input);
  virtual ~OutputManager() = default;
  void dump() const;
};

}  // namespace core
}  // namespace gryphon

#endif