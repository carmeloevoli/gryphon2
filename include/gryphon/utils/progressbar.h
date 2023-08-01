/* ProgressBar.h
 * Adapted from CRPropa3
 * https://crpropa.desy.de/
 */

#ifndef GRYPHON_UTIL_PROGRESSBAR_H
#define GRYPHON_UTIL_PROGRESSBAR_H

#include <ctime>
#include <memory>
#include <mutex>
#include <string>

namespace gryphon {
namespace utils {

class ProgressBar {
 private:
  unsigned long _steps;
  unsigned long _currentCount;
  unsigned long _maxbarLength;
  unsigned long _nextStep;
  unsigned long _updateSteps;
  time_t _startTime;
  std::string stringTmpl;
  std::string arrow;
  bool mutexSet;
  std::shared_ptr<std::mutex> _mutex;

 public:
  /// Initialize a ProgressBar with [steps] number of steps, updated at
  /// [updateSteps] intervalls
  ProgressBar(unsigned long steps = 0, unsigned long updateSteps = 100);
  /// Set mutex when used inside std::thread
  void setMutex(std::shared_ptr<std::mutex> mutex);
  /// Print a given title
  void start(const std::string &title);

  /// update the progressbar
  /// should be called steps times in a loop
  void update();

  // sets the position of the bar to a given value
  void setPosition(unsigned long position);

  /// Mark the progressbar with an error
  void setError();
};

}  // namespace utils
}  // namespace gryphon

#endif  // GRYPHON_UTIL_PROGRESSBAR_H