#ifndef GRYPHON_UTILS_LOGGING_H
#define GRYPHON_UTILS_LOGGING_H

#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Appenders/ConsoleAppender.h>
#include <plog/Appenders/RollingFileAppender.h>
#include <plog/Formatters/CsvFormatter.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Init.h>
#include <plog/Log.h>

#include "gryphon/utils/git_revision.h"

namespace gryphon {
namespace utils {

void startup_information();

}  // namespace utils
}  // namespace gryphon

#endif  // GRYPHON_UTILS_LOGGING_H