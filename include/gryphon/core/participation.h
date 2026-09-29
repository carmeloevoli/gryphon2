#ifndef GRYPHON_CORE_PARTICIPATION_H
#define GRYPHON_CORE_PARTICIPATION_H

#include <cmath>
#include <stdexcept>

namespace gryphon::core {

// Streaming participation number (sum I)^2 / sum I^2. Rescaling by the
// largest contribution avoids squaring tiny dimensional fluxes. No source
// catalogue or per-source spectrum needs to be retained by this accumulator.
class Participation {
 public:
  void add(double weight) {
    if (!std::isfinite(weight) || weight < 0.)
      throw std::invalid_argument("source weights must be finite and nonnegative");
    if (weight == 0.) return;
    if (weight > scale_) {
      const double ratio = scale_ / weight;
      sum_ *= ratio;
      squares_ *= ratio * ratio;
      scale_ = weight;
    }
    const double scaled = weight / scale_;
    sum_ += scaled;
    squares_ += scaled * scaled;
  }

  double effectiveCount() const { return squares_ > 0. ? sum_ * sum_ / squares_ : 0.; }
  double maxFraction() const { return sum_ > 0. ? 1. / sum_ : 0.; }

 private:
  double scale_ = 0.;
  double sum_ = 0.;
  double squares_ = 0.;
};

}  // namespace gryphon::core
#endif
