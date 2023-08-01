#ifndef GRYPHON_UTILS_GIT_REVISION_H
#define GRYPHON_UTILS_GIT_REVISION_H

#include <string>

std::string git_sha1();
std::string get_version();
bool git_has_local_changes();

#endif  // GRYPHON_UTILS_GIT_REVISION_H