#ifndef GRYPHON_UTILS_CLI_H
#define GRYPHON_UTILS_CLI_H

#include <optional>
#include <string>

#include "gryphon/core/input.h"

namespace gryphon {
namespace utils {

// Command line shared by every gryphon executable:
//
//     ./<program> <config.ini> [--seed N] [--outdir DIR]
//
// The configuration file carries the physics, the seed selects the
// realization, and the output directory is where runs are collected.
struct CommandLine {
  std::string config;
  std::optional<unsigned long int> seed;
  std::string outdir;
};

// Parses the command line, throwing std::invalid_argument with a usage string
// on anything unexpected. Returns std::nullopt when help was requested.
std::optional<CommandLine> parseCommandLine(int argc, char* argv[], const std::string& program);

// The input described by the command line: file, then seed override.
core::Input makeInput(const CommandLine& args);

}  // namespace utils
}  // namespace gryphon

#endif  // GRYPHON_UTILS_CLI_H
