#ifndef GRYPHON_GALAXY_FACTORY_H
#define GRYPHON_GALAXY_FACTORY_H

#include <memory>

#include "gryphon/core/input.h"
#include "gryphon/galaxy/galaxy.h"

namespace gryphon {
namespace galaxy {

std::shared_ptr<Galaxy> makeGalaxy(const core::Input& input);

}  // namespace galaxy
}  // namespace gryphon

#endif  // GRYPHON_GALAXY_FACTORY_H
