#include "gryphon/utils/cli.h"

#include <iostream>
#include <stdexcept>
#include <vector>

#include "gryphon/core/run.h"
#include "gryphon/utils/io.h"

namespace gryphon {
namespace utils {

namespace {

std::string usage(const std::string& program) {
  return "usage: ./" + program +
         " <config.ini> [--seed N] [--outdir DIR]\n"
         "  <config.ini>   parameter file; its name sets the simulation name\n"
         "  --seed N       overrides the seed in the parameter file\n"
         "  --outdir DIR   base output directory (default: " +
         std::string(core::RunOutput::kDefaultBaseDirectory) + ")";
}

[[noreturn]] void fail(const std::string& program, const std::string& message) {
  throw std::invalid_argument(message + "\n" + usage(program));
}

std::string requireValue(const std::string& program, const std::vector<std::string>& args,
                         size_t& index) {
  const std::string& flag = args[index];
  if (index + 1 >= args.size()) fail(program, "missing value for " + flag);
  return args[++index];
}

}  // namespace

std::optional<CommandLine> parseCommandLine(int argc, char* argv[], const std::string& program) {
  const std::vector<std::string> args(argv + 1, argv + argc);

  CommandLine parsed;
  parsed.outdir = core::RunOutput::kDefaultBaseDirectory;

  for (size_t i = 0; i < args.size(); ++i) {
    const std::string& arg = args[i];

    if (arg == "-h" || arg == "--help") {
      std::cout << usage(program) << "\n";
      return std::nullopt;
    } else if (arg == "--seed") {
      parsed.seed = parseSeed(requireValue(program, args, i).c_str());
    } else if (arg == "--outdir") {
      parsed.outdir = requireValue(program, args, i);
    } else if (arg == "--config") {
      parsed.config = requireValue(program, args, i);
    } else if (!arg.empty() && arg[0] == '-') {
      fail(program, "unknown option '" + arg + "'");
    } else if (parsed.config.empty()) {
      parsed.config = arg;
    } else {
      fail(program, "unexpected argument '" + arg + "'");
    }
  }

  if (parsed.config.empty()) fail(program, "no configuration file given");

  return parsed;
}

core::Input makeInput(const CommandLine& args) {
  core::Input input(args.config);
  if (args.seed) input.set_seed(*args.seed);
  return input;
}

}  // namespace utils
}  // namespace gryphon
