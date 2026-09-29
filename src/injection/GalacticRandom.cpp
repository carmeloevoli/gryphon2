#include "gryphon/injection/GalacticRandom.h"

#include <gsl/gsl_sf_gamma.h>

#include <stdexcept>

#include "gryphon/utils/numeric.h"

namespace gryphon {

namespace {

double pickSnEnergy(RandomNumberGenerator& rng) {
  auto logEnergy = rng.normal(0., 0.54);  // log10 scatter, in dex; median multiplier = 1
  return std::max(std::pow(10., logEnergy), 0.);
}

double pickSlope(double mean, double sdev, RandomNumberGenerator& rng, bool uncut) {
  auto slope = rng.normal(mean, sdev);
  if (uncut) {
    // Rejection sampling of the conditional Gaussian, not clipping at 2.
    // The parent mean is validated >2, so acceptance is at least one half.
    while (slope <= 2.) slope = rng.normal(mean, sdev);
    return slope;
  }
  return std::max(0.01, slope);
}

}  // namespace

namespace injection {

GalacticRandomSpectrum::GalacticRandomSpectrum(const core::Input& in, RandomNumberGenerator& rng)
    : InjectionSpectrum(in), m_E0(in.E_0()), m_Emin(10. * cgs::GeV) {
  m_Emax = in.injEmax();
  if (m_Emax == 0. && !(in.injSlope() > 2.)) {
    throw std::invalid_argument("Uncut GalacticRandom requires parent injSlope > 2");
  }
  m_alpha = (in.doVarySlope())
                ? pickSlope(in.injSlope(), in.injSlopeSigma(), rng, m_Emax == 0.)
                : in.injSlope();
  m_crenergy = in.injEfficiency() * cgs::E_SN * ((in.doVaryEnergy()) ? pickSnEnergy(rng) : 1.);
  m_Q0 = source_normalization();
}

double GalacticRandomSpectrum::source_normalization() const {
  const double prefactor = m_crenergy / pow2(m_E0);
  if (m_Emax == 0.) {
    // Exact Emax -> infinity limit of the energy integral above 10 GeV.
    // Zero is the explicit no-cutoff sentinel, not a very large finite cutoff.
    if (!(m_alpha > 2.)) {
      throw std::invalid_argument("Uncut GalacticRandom energy integral requires index > 2");
    }
    return prefactor * (m_alpha - 2.) * std::pow(m_Emin / m_E0, m_alpha - 2.);
  }
  if (m_Emax <= m_Emin) {
    throw std::invalid_argument(
        "GalacticRandomSpectrum requires injEmax = 0 or > 10 GeV for normalization");
  }

  // Integrate the same exponential cutoff used by get(), from 10 GeV to
  // infinity. Recur from a nonnegative gamma-function order: GSL 2.8's
  // small-x negative-order branch is inaccurate near the reference alpha.
  const double max_ratio = m_Emax / m_E0;
  const double exponent = 2. - m_alpha;
  const double x = m_Emin / m_Emax;
  const int shifts = (exponent < 0.) ? static_cast<int>(std::ceil(-exponent)) : 0;
  double gamma = gsl_sf_gamma_inc(exponent + shifts, x);
  for (int step = 1; step <= shifts; ++step) {
    const double order = exponent + (shifts - step);
    if (std::abs(order) < 1e-5) {
      // Downward recurrence loses precision close to order zero. These rare
      // near-integer indices use direct quadrature in log(E/E0) instead.
      const auto integrand = [&](double logRatio) {
        return std::exp(exponent * logRatio - std::exp(logRatio) / max_ratio);
      };
      const double integral = utils::QAGIntegration<double>(
          integrand, std::log(m_Emin / m_E0), std::log(100. * max_ratio), 1000, 1e-10);
      return prefactor / integral;
    }
    gamma = (gamma - std::exp(order * std::log(x) - x)) / order;
  }
  const double integral = std::pow(max_ratio, exponent) * gamma;
  return prefactor / integral;
}

double GalacticRandomSpectrum::get(double E) const {
  auto value = m_Q0 * std::pow(E / m_E0, -m_alpha);
  if (m_Emax > 0.) value *= std::exp(-(E / m_Emax));
  return value;
}

}  // namespace injection
}  // namespace gryphon
