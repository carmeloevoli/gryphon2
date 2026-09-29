#ifndef GRYPHON_CORE_RUN_H
#define GRYPHON_CORE_RUN_H

#include <filesystem>
#include <string>

#include "gryphon/core/input.h"
#include "gryphon/utils/io.h"

namespace gryphon {
namespace core {

// Output directory of one simulated model.
//
//     <base>/<simname>/params.ini            the resolved input, re-runnable
//     <base>/<simname>/<stem>_<seed>.txt     one table per seed
//
// Every table opens with a provenance header naming the code version, the
// parameters it came from, and its columns, so that a figure can always be
// traced back to the run that produced it.
class RunOutput {
 public:
  static constexpr const char* kDefaultBaseDirectory = "runs";

  RunOutput(const Input& input, const std::filesystem::path& baseDirectory);
  explicit RunOutput(const Input& input) : RunOutput(input, kDefaultBaseDirectory) {}

  const std::filesystem::path& directory() const noexcept { return m_directory; }

  // Full path of a table, e.g. open("flux") -> <dir>/flux_000003.txt
  std::filesystem::path path(const std::string& stem) const;

  // `columns` describes the table, one entry per column, separated by '|':
  //     "E [GeV] | I [GeV^-1 m^-2 s^-1 sr^-1]"
  utils::OutputFile open(const std::string& stem, const std::string& columns) const;

  std::string header(const std::string& columns) const;

 private:
  Input m_input;
  std::filesystem::path m_directory;
};

}  // namespace core
}  // namespace gryphon

#endif  // GRYPHON_CORE_RUN_H
