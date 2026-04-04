#ifndef GRYPHON_INJECTION_INJECTIONSPECTRUM_H
#define GRYPHON_INJECTION_INJECTIONSPECTRUM_H

#include "gryphon/core/input.h"

namespace gryphon {
namespace injection {

class InjectionSpectrum {
 public:
  explicit InjectionSpectrum(const core::Input& in) : m_in(in) {}
  virtual ~InjectionSpectrum() = default;

  virtual double get(double E) const = 0;

 protected:
  const core::Input& m_in;
};

}  // namespace injection
}  // namespace gryphon

#endif  // GRYPHON_INJECTION_INJECTIONSPECTRUM_H
