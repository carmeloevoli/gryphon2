#include "gryphon.h"

using namespace gryphon;

int main() {
  try {
    utils::startup_information();
  } catch (std::exception& e) {
    LOGE << "!Fatal Error: " << e.what();
  }
  return 0;
}