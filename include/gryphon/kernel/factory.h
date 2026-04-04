#ifndef GRYPHON_KERNEL_FACTORY_H
#define GRYPHON_KERNEL_FACTORY_H

#include <memory>

#include "gryphon/core/input.h"
#include "gryphon/kernel/greenkernel.h"

namespace gryphon {
namespace kernel {

std::shared_ptr<const GreenKernel> makeGreenKernel(const core::Input& in);

}  // namespace kernel
}  // namespace gryphon

#endif  // GRYPHON_KERNEL_FACTORY_H
