#include "gryphon/galaxy/factory.h"

#include <stdexcept>

#include "gryphon/galaxy/Steiman2010.h"
#include "gryphon/galaxy/jelly.h"
#include "gryphon/galaxy/uniform.h"

namespace gryphon {
namespace galaxy {

std::shared_ptr<Galaxy> makeGalaxy(const core::Input& input) {
  input.validate();
  switch (input.spiralModel()) {
    case SpiralModel::Uniform:
      return std::make_shared<GalaxyUniform>(input);
    case SpiralModel::Jelly:
      return std::make_shared<GalaxyJelly>(input);
    case SpiralModel::Steiman2010:
      return std::make_shared<GalaxySteiman2010>(input);
    default:
      throw std::invalid_argument("Spiral model not implemented yet");
  }
}

}  // namespace galaxy
}  // namespace gryphon
