#ifndef GRYPHON_UTILS_IO_H
#define GRYPHON_UTILS_IO_H

#include <fstream>
#include <string>
#include <vector>

namespace gryphon {
namespace utils {

// String Utilities
std::string removeExtensionIniFilename(std::string inputFilename);

// Files
size_t countFileLines(const std::string& filename);
bool fileExists(const std::string& filename);
std::vector<std::string> split(std::string s, std::string delimiter = " ");
std::vector<double> loadRow(std::string filePath, size_t iRow, std::string delimiter = " ");
std::vector<std::vector<double> > loadFileByRow(std::string filePath, std::string delimiter = " ");

// Output file
class OutputFile {
  std::string filename;
  std::ofstream out;

 public:
  OutputFile(const std::string& name);
  ~OutputFile();

  template <typename T>
  OutputFile& operator<<(const T& value) {
    out << value;
    return *this;
  }
};

}  // namespace utils
}  // namespace gryphon

#endif  // GRYPHON_UTILS_IO_H