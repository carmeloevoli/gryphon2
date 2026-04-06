#ifndef GRYPHON_INJECTION_FACTORY_H
#define GRYPHON_INJECTION_FACTORY_H

#include "gryphon/core/event.h"
#include "gryphon/core/input.h"
#include "gryphon/injection/InjectionSpectrum.h"
#include "gryphon/utils/random.h"

namespace gryphon {
namespace injection {

InjectionSpectrumPtr makeInjectionSpectrum(const core::Input& in, RandomNumberGenerator& rng);
InjectionSpectra makeInjectionSpectra(const core::Input& in, const core::Events& events,
                                      RandomNumberGenerator& rng);

}  // namespace injection
}  // namespace gryphon

#endif  // GRYPHON_INJECTION_FACTORY_H
