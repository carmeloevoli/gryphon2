#include "gryphon/utils/numeric.h"

#include <algorithm>
#include <cmath>

#include "gryphon/utils/spline.h"

namespace gryphon {
namespace utils {

double halo_function(double l2, double H, double z, double zs, double rel_error) {
  constexpr double kPi = 3.141592653589793238462643383279502884;
  if (!(l2 > 0.) || !(H > 0.)) return 0.;
  const double tol = std::max(rel_error, 1e-15);
  const double lambdaOverH = std::sqrt(l2) / H;

  auto G1D = [](double x, double s2) { return std::exp(-(x * x) / s2); };

  // For large diffusion lengths, the image expansion suffers from strong
  // cancellation; the eigenmode expansion converges much more stably.
  if (lambdaOverH > 3.) {
    const size_t n_max = 2000;
    const double pref = std::sqrt(kPi * l2) / (2. * H);
    const double x = (z + H) / (2. * H);
    const double y = (zs + H) / (2. * H);
    const double decay = kPi * kPi * l2 / (16. * H * H);

    double sum = 0.;
    double c = 0.;
    for (size_t n = 1; n <= n_max; ++n) {
      const double nn = static_cast<double>(n);
      const double term =
          std::exp(-nn * nn * decay) * std::sin(nn * kPi * x) * std::sin(nn * kPi * y);

      // Kahan compensated summation for cancellation-prone alternating tails.
      const double yk = term - c;
      const double tk = sum + yk;
      c = (tk - sum) - yk;
      sum = tk;

      if (std::abs(term) <= tol * std::max(1., std::abs(sum))) break;
    }
    return std::max(pref * sum, 0.);
  }

  // Image sum is efficient and accurate for small-to-intermediate diffusion
  // lengths.
  double f = G1D(z - zs, l2);
  double c = 0.;
  const size_t n_max = 2000;
  for (size_t n = 1; n <= n_max; ++n) {
    const double sign = (n % 2 == 0) ? 1. : -1.;
    const double nn = static_cast<double>(n);
    const double zn_plus = sign * zs + 2. * nn * H;
    const double zn_minus = sign * zs - 2. * nn * H;
    const double term = sign * (G1D(z - zn_plus, l2) + G1D(z - zn_minus, l2));

    const double yk = term - c;
    const double tk = f + yk;
    c = (tk - f) - yk;
    f = tk;

    if (std::abs(term) <= tol * std::max(1., std::abs(f))) break;
  }

  return std::max(f, 0.);
}

#define index(i, j) ((j) + (i) * Y.size())

double interpolate(double x, const std::vector<double> &X, const std::vector<double> &Y) {
  std::vector<double>::const_iterator it = std::upper_bound(X.begin(), X.end(), x);
  if (it == X.begin()) return 0;
  if (it == X.end()) return 0;

  const size_t i = it - X.begin() - 1;
  return Y[i] + (x - X[i]) * (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]);
}

double interpolateEquidistant(double x, double lo, double hi, const std::vector<double> &Y) {
  if (x <= lo) return 0;
  if (x >= hi) return 0;

  const double dx = (hi - lo) / (Y.size() - 1);
  const double p = (x - lo) / dx;
  const size_t i = std::floor(p);
  return Y[i] + (p - i) * (Y[i + 1] - Y[i]);
}

double cspline(double x, const std::vector<double> &X, const std::vector<double> &Y) {
  if (x < X.front()) return 0;
  if (x > X.back()) return 0;
  tk::spline s(X, Y, tk::spline::cspline);
  return s(x);
}

double interpolate2d(double x, double y, const std::vector<double> &X, const std::vector<double> &Y,
                     const std::vector<double> &Z) {
  std::vector<double>::const_iterator itx = std::upper_bound(X.begin(), X.end(), x);
  std::vector<double>::const_iterator ity = std::upper_bound(Y.begin(), Y.end(), y);

  if (x > X.back() || x < X.front()) return 0;
  if (y > Y.back() || y < Y.front()) return 0;

  if (itx == X.begin() && ity == Y.begin()) return Z.front();
  if (itx == X.end() && ity == Y.end()) return Z.back();

  const size_t i = itx - X.begin() - 1;
  const size_t j = ity - Y.begin() - 1;

  const double Q11 = Z[index(i, j)];
  const double Q12 = Z[index(i, j + 1)];
  const double Q21 = Z[index(i + 1, j)];
  const double Q22 = Z[index(i + 1, j + 1)];

  const double R1 =
      ((X[i + 1] - x) / (X[i + 1] - X[i])) * Q11 + ((x - X[i]) / (X[i + 1] - X[i])) * Q21;
  const double R2 =
      ((X[i + 1] - x) / (X[i + 1] - X[i])) * Q12 + ((x - X[i]) / (X[i + 1] - X[i])) * Q22;

  return ((Y[j + 1] - y) / (Y[j + 1] - Y[j])) * R1 + ((y - Y[j]) / (Y[j + 1] - Y[j])) * R2;
}

}  // namespace utils
}  // namespace gryphon
