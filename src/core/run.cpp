#include "gryphon/core/run.h"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>

#include "gryphon/utils/git_revision.h"
#include "gryphon/utils/logging.h"

namespace gryphon {
namespace core {

namespace fs = std::filesystem;

namespace {

constexpr int kSeedDigits = 6;

std::string timestamp() {
  const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::tm parts{};
#ifdef _WIN32
  gmtime_s(&parts, &now);
#else
  gmtime_r(&now, &parts);
#endif
  std::ostringstream stream;
  stream << std::put_time(&parts, "%Y-%m-%dT%H:%M:%SZ");
  return stream.str();
}

std::string paddedSeed(unsigned long int seed) {
  std::ostringstream stream;
  stream << std::setfill('0') << std::setw(kSeedDigits) << seed;
  return stream.str();
}

}  // namespace

RunOutput::RunOutput(const Input& input, const fs::path& baseDirectory)
    : m_input(input), m_directory(baseDirectory / input.simname()) {
  m_input.validate();

  std::error_code error;
  fs::create_directories(m_directory, error);
  if (error && !fs::is_directory(m_directory)) {
    throw std::runtime_error("could not create run directory '" + m_directory.string() +
                             "': " + error.message());
  }

  m_input.write_params_file((m_directory / "params.ini").string());
  LOGD << "run directory : " << m_directory.string();
}

fs::path RunOutput::path(const std::string& stem) const {
  return m_directory / (stem + "_" + paddedSeed(m_input.seed()) + ".txt");
}

std::string RunOutput::header(const std::string& columns) const {
  std::ostringstream out;
  out << "# gryphon = " << get_version() << "\n";
  out << "# git = " << git_sha1() << (git_has_local_changes() ? " (dirty)" : "") << "\n";
  out << "# date = " << timestamp() << "\n";
  out << "# simname = " << m_input.simname() << "\n";
  out << "# seed = " << m_input.seed() << "\n";
  if (!m_input.configFile().empty()) out << "# config = " << m_input.configFile() << "\n";
  out << "# params = params.ini\n";
  out << "# columns: " << columns << "\n";
  return out.str();
}

utils::OutputFile RunOutput::open(const std::string& stem, const std::string& columns) const {
  return utils::OutputFile(path(stem), header(columns));
}

}  // namespace core
}  // namespace gryphon
