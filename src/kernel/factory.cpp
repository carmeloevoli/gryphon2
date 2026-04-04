#include "gryphon/kernel/factory.h"

#include <stdexcept>

#include "gryphon/kernel/DiffusionLossesKernel.h"
#include "gryphon/kernel/PureDiffusionKernel.h"

namespace gryphon {
namespace kernel {

std::shared_ptr<const GreenKernel> makeGreenKernel(const core::Input& in) {
  in.validate();
  switch (in.transportModel()) {
    case TransportModel::PureDiffusion:
      return std::make_shared<PureDiffusionKernel>(in);
    case TransportModel::DiffusionLosses:
      return std::make_shared<DiffusionLossesKernel>(in);
    default:
      throw std::invalid_argument("Transport model not implemented yet");
  }
}

}  // namespace kernel
}  // namespace gryphon
