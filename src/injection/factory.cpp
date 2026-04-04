#include "gryphon/injection/factory.h"

#include <stdexcept>

#include "gryphon/injection/GalacticRandom.h"
#include "gryphon/injection/SinglePowerLaw.h"

namespace gryphon {
namespace injection {

std::shared_ptr<InjectionSpectrum> makeInjectionSpectrum(const core::Input& in,
                                                         RandomNumberGenerator& rng) {
  in.validate();
  switch (in.injectionModel()) {
    case InjectionModel::SinglePowerLaw:
      return std::make_shared<SinglePowerLawSpectrum>(in);
    case InjectionModel::GalacticRandom:
      return std::make_shared<GalacticRandomSpectrum>(in, rng);
    default:
      throw std::invalid_argument("Injection model not implemented yet");
  }
}

}  // namespace injection
}  // namespace gryphon
