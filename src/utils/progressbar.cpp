#include "gryphon/utils/progressbar.h"

#include <cstdio>
#include <iostream>
#include <utility>

#include "gryphon/utils/logging.h"

namespace gryphon {
namespace utils {

/// Initialize a ProgressBar with [steps] number of steps, updated at
/// [updateSteps] intervalls
ProgressBar::ProgressBar(unsigned long steps, unsigned long updateSteps)
    : _steps(steps),
      _currentCount(0),
      _maxbarLength(10),
      _nextStep(1),
      _updateSteps(updateSteps),
      _startTime(0),
      mutexSet(false) {
  if (_updateSteps > _steps) _updateSteps = _steps;
  arrow.append(">");
}

void ProgressBar::setMutex(std::shared_ptr<std::mutex> mutex) {
  _mutex = std::move(mutex);
  mutexSet = true;
}

void ProgressBar::start(const std::string &title) {
  _startTime = time(NULL);
  std::string s = ctime(&_startTime);
  s.erase(s.end() - 1, s.end());
  stringTmpl = "  Started ";
  stringTmpl.append(s);
  stringTmpl.append(" : [%-10s] %3i%%    %s: %02i:%02i:%02i %s\r");
  LOGD << title;
}
/// update the progressbar
/// should be called steps times in a loop
void ProgressBar::update() {
  if (mutexSet) {
    std::lock_guard<std::mutex> guard(*_mutex);
    _currentCount++;
    if (_currentCount == _nextStep || _currentCount == _steps || _currentCount == 1000) {
      _nextStep += long(_steps / float(_updateSteps));
      setPosition(_currentCount);
    }
  } else {
    _currentCount++;
    if (_currentCount == _nextStep || _currentCount == _steps || _currentCount == 1000) {
      _nextStep += long(_steps / float(_updateSteps));
      setPosition(_currentCount);
    }
  }
}

void ProgressBar::setPosition(unsigned long position) {
  int percentage = int(100 * (position / float(_steps)));
  time_t currentTime = time(NULL);
  if (position < _steps) {
    if (arrow.size() <= (_maxbarLength) * (position) / (_steps)) arrow.insert(0, "=");
    time_t tElapsed = currentTime - _startTime;
    float tToGo = (_steps - position) * tElapsed / position;
    std::printf(stringTmpl.c_str(), arrow.c_str(), percentage, "Finish in", int(tToGo / 3600),
                (int(tToGo) % 3600) / 60, int(tToGo) % 60, "");
    fflush(stdout);
  } else {
    float tElapsed = currentTime - _startTime;
    std::string s = " - Finished at ";
    s.append(ctime(&currentTime));
    char fs[255];
    std::snprintf(fs, 100, "%c[%d;%dm Finished %c[%dm", 27, 1, 32, 27, 0);
    std::printf(stringTmpl.c_str(), fs, 100, "Needed", int(tElapsed / 3600),
                (int(tElapsed) % 3600) / 60, int(tElapsed) % 60, s.c_str());
  }
}

/// Mark the progressbar with an error
void ProgressBar::setError() {
  time_t currentTime = time(NULL);
  _currentCount++;
  time_t tElapsed = currentTime - _startTime;
  std::string s = " - Finished at ";
  s.append(ctime(&currentTime));
  char fs[255];
  std::snprintf(fs, 100, "%c[%d;%dm  ERROR   %c[%dm", 27, 1, 31, 27, 0);
  std::printf(stringTmpl.c_str(), fs, _currentCount, "Needed", int(tElapsed / 3600),
              (int(tElapsed) % 3600) / 60, int(tElapsed) % 60, s.c_str());
}

}  // namespace utils
}  // namespace gryphon