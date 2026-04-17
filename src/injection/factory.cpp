#include "gryphon/injection/factory.h"

#include <stdexcept>

#include "gryphon/injection/GalacticRandom.h"
#include "gryphon/injection/PWN.h"
#include "gryphon/injection/RandomEmax.h"
#include "gryphon/injection/SinglePowerLaw.h"
#include "gryphon/injection/YoungPulsars.h"

namespace gryphon {
namespace injection {

InjectionSpectrumPtr makeInjectionSpectrum(const core::Input& in, RandomNumberGenerator& rng) {
  in.validate();
  switch (in.injectionModel()) {
    case InjectionModel::SinglePowerLaw:
      return std::make_shared<SinglePowerLawSpectrum>(in);
    case InjectionModel::GalacticRandom:
      return std::make_shared<GalacticRandomSpectrum>(in, rng);
    case InjectionModel::RandomEmax:
      return std::make_shared<RandomEmaxSpectrum>(in, rng);
    case InjectionModel::PWN:
      return std::make_shared<PWNSpectrum>(in, rng);
    case InjectionModel::YoungPulsars:
      return std::make_shared<YoungPulsarsSpectrum>(in, rng);
    default:
      throw std::invalid_argument("Injection model not implemented yet");
  }
}

InjectionSpectra makeInjectionSpectra(const core::Input& in, const core::Events& events,
                                      RandomNumberGenerator& rng) {
  in.validate();

  InjectionSpectra spectra;
  spectra.reserve(events.size());

  switch (in.injectionModel()) {
    case InjectionModel::SinglePowerLaw: {
      const auto shared_spectrum = std::make_shared<SinglePowerLawSpectrum>(in);
      for (const auto& event : events) {
        spectra.push_back(event ? shared_spectrum : nullptr);
      }
      return spectra;
    }
    case InjectionModel::GalacticRandom: {
      if (!in.doVaryEnergy() && !in.doVarySlope()) {
        const auto shared_spectrum = std::make_shared<GalacticRandomSpectrum>(in, rng);
        for (const auto& event : events) {
          spectra.push_back(event ? shared_spectrum : nullptr);
        }
        return spectra;
      }

      for (const auto& event : events) {
        if (!event) {
          spectra.push_back(nullptr);
          continue;
        }
        spectra.push_back(std::make_shared<GalacticRandomSpectrum>(in, rng));
      }
      return spectra;
    }
    case InjectionModel::RandomEmax: {
      if (in.injEmaxSigmaDex() <= 0.) {
        const auto shared_spectrum = std::make_shared<RandomEmaxSpectrum>(in, rng);
        for (const auto& event : events) {
          spectra.push_back(event ? shared_spectrum : nullptr);
        }
        return spectra;
      }

      for (const auto& event : events) {
        if (!event) {
          spectra.push_back(nullptr);
          continue;
        }
        spectra.push_back(std::make_shared<RandomEmaxSpectrum>(in, rng));
      }
      return spectra;
    }
    case InjectionModel::PWN: {
      for (const auto& event : events) {
        if (!event) {
          spectra.push_back(nullptr);
          continue;
        }
        spectra.push_back(std::make_shared<PWNSpectrum>(in, rng));
      }
      return spectra;
    }
    case InjectionModel::YoungPulsars: {
      for (const auto& event : events) {
        if (!event) {
          spectra.push_back(nullptr);
          continue;
        }
        spectra.push_back(std::make_shared<YoungPulsarsSpectrum>(in, rng));
      }
      return spectra;
    }
    default:
      throw std::invalid_argument("Injection model not implemented yet");
  }
}

}  // namespace injection
}  // namespace gryphon
