#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

#include "gryphon.h"
#include "gtest/gtest.h"

namespace gryphon {
namespace {

namespace fs = std::filesystem;

std::string makeTemporaryDirectory() {
  char pattern[] = "/tmp/gryphon_io_XXXXXX";
  char* directory = ::mkdtemp(pattern);
  if (directory == nullptr) throw std::runtime_error("failed to create temporary directory");
  return std::string(directory);
}

class ScopedWorkingDirectory {
 public:
  explicit ScopedWorkingDirectory(const fs::path& path) : previous_(fs::current_path()) {
    fs::current_path(path);
  }

  ~ScopedWorkingDirectory() {
    std::error_code error;
    fs::current_path(previous_, error);
  }

 private:
  fs::path previous_;
};

}  // namespace

TEST(OutputFile, CreatesOutputDirectoryAndWritesContent) {
  const fs::path tempDirectory = makeTemporaryDirectory();
  {
    ScopedWorkingDirectory cwd(tempDirectory);

    {
      utils::OutputFile out("sample.txt");
      out << "alpha\n";
      out << 42 << "\n";
    }

    ASSERT_TRUE(utils::fileExists("output/sample.txt"));
    ASSERT_TRUE(fs::exists("output/sample.txt"));

    std::ifstream file("output/sample.txt");
    ASSERT_TRUE(file.is_open());

    std::string firstLine;
    std::string secondLine;
    ASSERT_TRUE(static_cast<bool>(std::getline(file, firstLine)));
    ASSERT_TRUE(static_cast<bool>(std::getline(file, secondLine)));
    EXPECT_EQ(firstLine, "alpha");
    EXPECT_EQ(secondLine, "42");

    EXPECT_TRUE(fs::remove("output/sample.txt"));
    EXPECT_TRUE(fs::remove("output"));
  }
  EXPECT_TRUE(fs::remove(tempDirectory));
}

TEST(ParseSeed, AcceptsValidUnsignedIntegers) {
  EXPECT_EQ(utils::parseSeed("0"), 0u);
  EXPECT_EQ(utils::parseSeed("42"), 42u);
}

TEST(ParseSeed, RejectsInvalidValues) {
  EXPECT_THROW(utils::parseSeed(nullptr), std::invalid_argument);
  EXPECT_THROW(utils::parseSeed(""), std::invalid_argument);
  EXPECT_THROW(utils::parseSeed("-1"), std::invalid_argument);
  EXPECT_THROW(utils::parseSeed("abc"), std::invalid_argument);
  EXPECT_THROW(utils::parseSeed("12abc"), std::invalid_argument);
}

TEST(FileLoading, IgnoresBlankLinesAndComments) {
  const fs::path tempDirectory = makeTemporaryDirectory();
  {
    ScopedWorkingDirectory cwd(tempDirectory);

    std::ofstream file("rows.csv");
    ASSERT_TRUE(file.is_open());
    file << "# comment\n";
    file << "\n";
    file << "1,2,3\n";
    file << "   \n";
    file << "4,5,6\n";
    file.close();

    EXPECT_EQ(utils::countFileLines("rows.csv"), 2u);

    const auto row = utils::loadRow("rows.csv", 1, ",");
    ASSERT_EQ(row.size(), 3u);
    EXPECT_DOUBLE_EQ(row[0], 4.0);
    EXPECT_DOUBLE_EQ(row[1], 5.0);
    EXPECT_DOUBLE_EQ(row[2], 6.0);

    const auto rows = utils::loadFileByRow("rows.csv", ",");
    ASSERT_EQ(rows.size(), 2u);
    EXPECT_DOUBLE_EQ(rows[0][0], 1.0);
    EXPECT_DOUBLE_EQ(rows[1][2], 6.0);

    EXPECT_TRUE(fs::remove("rows.csv"));
  }
  EXPECT_TRUE(fs::remove(tempDirectory));
}

TEST(XsecsLoading, MissingFileThrows) {
  EXPECT_THROW(Orusa2022::XSECS("does-not-exist.txt"), std::runtime_error);
}

TEST(XsecsLoading, BundledTableParsesSuccessfully) {
  const Orusa2022::PositronXsecs xsecs;
  const double value = xsecs.get(Orusa2022::Channel::pp, 10.0 * cgs::GeV, 1.0 * cgs::GeV);

  EXPECT_TRUE(std::isfinite(value));
  EXPECT_GE(value, 0.0);
}

}  // namespace gryphon
