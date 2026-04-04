#ifndef GRYPHON_INJECTION_FACTORY_H
#define GRYPHON_INJECTION_FACTORY_H

#include <memory>

#include "gryphon/core/input.h"
#include "gryphon/injection/InjectionSpectrum.h"
#include "gryphon/utils/random.h"

namespace gryphon {
namespace injection {

std::shared_ptr<InjectionSpectrum> makeInjectionSpectrum(const core::Input& in,
                                                         RandomNumberGenerator& rng);

}  // namespace injection
}  // namespace gryphon

#endif  // GRYPHON_INJECTION_FACTORY_H
